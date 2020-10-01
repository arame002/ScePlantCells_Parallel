
#include "SignalTissue.hpp"
#include <iomanip>
#include <sstream>

void SignalTissue::Cal_AllCellCenters()
{
    for (unsigned int i = 0 ; i< cells.size(); i++)
    {
        cells.at(i).Cal_Centroid() ;
    }
    
}
//---------------------------------------------------------------------------------------------
void SignalTissue::Cal_TissueCenter()
{
    double tissueCntrX = 0.0 ;
    double tissueCntrY = 0.0 ;
    for (unsigned int i=0; i< cells.size(); i++)
    {
        tissueCntrX += cells.at(i).centroid.at(0) ;
        tissueCntrY += cells.at(i).centroid.at(1) ;
    }
    tissueCenter.push_back(tissueCntrX/cells.size() ) ;
    tissueCenter.push_back(tissueCntrY/cells.size() ) ;
}
//---------------------------------------------------------------------------------------------
void SignalTissue::Find_AllCellNeighborCandidates ()
{
    for (unsigned int i=0; i < cells.size(); i++)
    {
        cells.at(i).Find_NghbrCandidate() ;
    }
}
//---------------------------------------------------------------------------------------------

void SignalTissue::Cal_AllCellCntrToCntr()
{
    for (unsigned int i=0 ; i< cells.size() ; i++ )
    {
        for (unsigned int j=0 ; j< cells.size(); j++)
        {
            cells.at(i).cntrToCntr.push_back(Dist2D(cells.at(i).centroid.at(0),cells.at(i).centroid.at(1),
                                                    cells.at(j).centroid.at(0),cells.at(j).centroid.at(1))) ;
            
        }
    }
}

//---------------------------------------------------------------------------------------------
void SignalTissue::FindInterfaceWithNeighbor()
{
    for (unsigned int i=0; i < cells.size(); i++)
    {
        for (unsigned int j= 0; j < cells.at(i).neighbors.size(); j++)
        {
            int cellIDNeighbor = cells.at(i).neighbors.at(j).CellID_Neighbor ;
            int neighborIDCell = cells.at(i).neighbors.at(j).NeighborID_Cell ;
            
            if (cells.at(i).neighbors.at(j).common_edge_ID == false )
            {
        
                for (unsigned int k =0; k < cells.at(i).neighbors.at(j).indexCellNode.size() ; k++)
                {
                    int nodeIDCell = cells.at(i).neighbors.at(j).indexCellNode.at(k) ;
                    int nodeIDNeighbor = cells.at(i).neighbors.at(j).indexMinDistNode.at(nodeIDCell) ;
                    double newX = 0.5*(cells.at(i).nodesX.at(nodeIDCell) + cells.at(cellIDNeighbor).nodesX.at(nodeIDNeighbor) ) ;
                    double newY = 0.5*(cells.at(i).nodesY.at(nodeIDCell) + cells.at(cellIDNeighbor).nodesY.at(nodeIDNeighbor) ) ;
                
                    cells.at(i).nodesXNew.at(nodeIDCell)= newX ;
                    cells.at(i).nodesYNew.at(nodeIDCell)= newY ;
                    cells.at(cellIDNeighbor).nodesXNew.at(nodeIDNeighbor) = newX ;
                    cells.at(cellIDNeighbor).nodesYNew.at(nodeIDNeighbor) = newY ;
        // we should work with either nodesXNew or intfX
                    
                    cells.at(i).neighbors.at(j).intfX.push_back(newX) ;
                    cells.at(i).neighbors.at(j).intfY.push_back(newY) ;
                    cells.at(cellIDNeighbor).neighbors.at(neighborIDCell).intfX.push_back(newX) ;
                    cells.at(cellIDNeighbor).neighbors.at(neighborIDCell).intfY.push_back(newY) ;
                
                  
                }
                cells.at(i).neighbors.at(j).common_edge_ID = true ;
                cells.at(cellIDNeighbor).neighbors.at(neighborIDCell).common_edge_ID = true ;
                
            }
        }
    }
    
}

//---------------------------------------------------------------------------------------------

vector<SignalCell> SignalTissue::ReadFile ( )
{
    double a,b,c ;
    double d ;
    vector<SignalCell> tmpCell ;
    string number = to_string(frameIndex) ;
    ifstream nodeData ("Locations_"+ number +".txt") ;
    if (nodeData.is_open())
    {
        cout << "file is open"<<endl ;
    }
    else{
        cout << "Error in file"<<endl ;
    }
    
    while (nodeData >> a >> b >> c >> d)
    {
        if (a==b)
        {
            SignalCell cell ;
            tmpCell.push_back(cell) ;
            tmpCell.back().cellID = a ;
            continue ;
        }
        
        if ( d == 1 )
        {
            
            tmpCell.back().nodesX.push_back(a) ;
            tmpCell.back().nodesY.push_back(b) ;
            tmpCell.back().nodesXNew.push_back(a) ;
            tmpCell.back().nodesYNew.push_back(b) ;
            tmpCell.back().noNeighboringNodesX.push_back(a) ;
            tmpCell.back().noNeighboringNodesY.push_back(b) ;
            
            
            //    cells.back().nodes.back().z = c ;
            //    cells.back().nodes.back().nodeType = d ;
        }
    }
    return tmpCell ;
}
//---------------------------------------------------------------------------------------------
vector<SignalCell> SignalTissue::ReadFile2 ( )
{
    int cellLayer = 0 ;
    int a,b ;
    double c,d ;
    double e,f ;
    vector<SignalCell> tmpCell ;
    string number = to_string(frameIndex) ;
    ifstream nodeData ("Locations_"+ number +".txt") ;
    if (nodeData.is_open())
    {
        cout << "file is open"<<endl ;
    }
    else{
        cout << "Error in file"<<endl ;
    }
    
    while (nodeData >> b >> a >> c >> d >> e >> f )
    {
        if (static_cast<int>( tmpCell.size() ) == b)
        {
            cellLayer = a ;
            SignalCell cell ;
            tmpCell.push_back(cell) ;
        }
        
        tmpCell.back().nodesX.push_back(c) ;
        tmpCell.back().nodesY.push_back(d) ;
        tmpCell.back().nodesXNew.push_back(c) ;
        tmpCell.back().nodesYNew.push_back(d) ;
        tmpCell.back().noNeighboringNodesX.push_back(c) ;
        tmpCell.back().noNeighboringNodesY.push_back(d) ;
        tmpCell.back().layer = cellLayer ;
            
            //    cells.back().nodes.back().z = c ;
            //    cells.back().nodes.back().nodeType = d ;
    }
    return tmpCell ;
}
//---------------------------------------------------------------------------------------------
vector<SignalCell> SignalTissue::ReadFile3 ( )
{
    int b ;
    double c,d ;
    vector<SignalCell> tmpCell ;
    string number = to_string(frameIndex) ;
    ifstream nodeData ("ExportCellProp_"+ number +".txt") ;
    if (nodeData.is_open())
    {
     //   cout << "file is open"<<endl ;
    }
    else{
        cout << "Error in file"<<endl ;
    }
    
    while ( nodeData  >> b >> c >> d )
    {
        if (static_cast<int>( tmpCell.end() - tmpCell.begin() ) == b )
        {
            SignalCell cell ;
            tmpCell.push_back(cell) ;
        }
        
        tmpCell.back().nodesX.push_back(c) ;
        tmpCell.back().nodesY.push_back(d) ;
        tmpCell.back().nodesXNew.push_back(c) ;
        tmpCell.back().nodesYNew.push_back(d) ;
        tmpCell.back().noNeighboringNodesX.push_back(c) ;
        tmpCell.back().noNeighboringNodesY.push_back(d) ;
        
        //    cells.back().nodes.back().z = c ;
        //    cells.back().nodes.back().nodeType = d ;
    }
    return tmpCell ;
}
//---------------------------------------------------------------------------------------------

void SignalTissue::Coupling (vector< vector<double> > locX , vector< vector<double> > locY)
{
    int cellSize = static_cast<int>( locX.size() ) ;
    for (unsigned int i =0 ; i< cellSize ; i++)
    {
        SignalCell tmpCell ;
        tmpCell.nodesX = locX.at(i) ;
        tmpCell.nodesXNew = locX.at(i) ;
        tmpCell.noNeighboringNodesX = locX.at(i) ;
        
        tmpCell.nodesY = locY.at(i) ;
        tmpCell.nodesYNew = locY.at(i) ;
        tmpCell.noNeighboringNodesY = locY.at(i) ;
        
        cells.push_back( tmpCell ) ;
    }
 
}

//---------------------------------------------------------------------------------------------

void SignalTissue::Find_AllCellNeighbors()
{
    for (unsigned int i=0; i < cells.size(); i++)
    {
        for (unsigned int j= 0; j < cells.at(i).nghbrCandidate.size(); j++)
        {
            int idNeighbor=cells.at(i).nghbrCandidate.at(j) ;
            
            cells.at(i).nodeDistoNghbrsCandidate.push_back (cells.at(i).Cal_NodeToNodeDist(cells.at(idNeighbor).nodesX ,
                                                                                           cells.at(idNeighbor).nodesY ) );
            
        }
        cells.at(i).Find_NghbrProperties() ;
    }
    Find_AllCell_NeighborID_Cell() ;
}

//---------------------------------------------------------------------------------------------
void SignalTissue::Find_AllCell_NeighborID_Cell()
{
    for (unsigned int i=0 ; i< cells.size(); i++)
    {
        for (unsigned int j=0; j< cells.at(i).neighbors.size(); j++)
        {
            int neighboringCell = cells.at(i).neighbors.at(j).CellID_Neighbor ;
            // search on all the neighbors of the neighboring cell
            for (unsigned int k=0; k < cells.at(neighboringCell).neighbors.size(); k++)
            {
              if (  cells.at(neighboringCell).neighbors.at(k).CellID_Neighbor == i )
              {
                  cells.at(i).neighbors.at(j).NeighborID_Cell = k ;
                  break ;
              }
            }
        }
    }
}
//---------------------------------------------------------------------------------------------

void SignalTissue::Cal_AllCellNewEdge()
{
    for (unsigned int i=0 ; i < cells.size(); i++)
    {
        cells.at(i).NewEdge() ;
    }
    
}
//---------------------------------------------------------------------------------------------

void SignalTissue::Find_CommonNeighbors()
{
    for (unsigned int i = 0; i < cells.size(); i++)
    {
        for (unsigned int j = 0; j < cells.at(i).neighbors.size(); j++)
        {
            int cellIDNeighbor = cells.at(i).neighbors.at(j).CellID_Neighbor ;
            int neighborIDCell = cells.at(i).neighbors.at(j).NeighborID_Cell ;
            
            if (i < cellIDNeighbor)
            {
                for (unsigned int k=j+1; k < cells.at(i).neighbors.size() ; k++)
                {
                    for (unsigned int l= 0; l < cells.at(cellIDNeighbor).neighbors.size() ; l++)
                    {
                        if ( cells.at(i).neighbors.at(k).CellID_Neighbor == cells.at(cellIDNeighbor).neighbors.at(l).CellID_Neighbor )
                        {
                            int cellID_Other = cells.at(i).neighbors.at(k).CellID_Neighbor ;
                            int nghbrIDCell_Other = cells.at(i).neighbors.at(k).NeighborID_Cell ;
                            int nghbrIDNeighboringCell_Other = cells.at(cellIDNeighbor).neighbors.at(l).NeighborID_Cell ;
                            
                            CommonNeighbor tmpAB ;
                            tmpAB.cellIDCommonNeighbor = cellID_Other ;
                            tmpAB.nghbrIDCell_CommonNeighbor = nghbrIDCell_Other ;
                            tmpAB.nghbrIDNeighboringCell_CommonNeighbor = nghbrIDNeighboringCell_Other ;
                            tmpAB.nghbrIDCommonNeighbor_Cell = k ;
                            tmpAB.nghbrIDCommonNeighbor_NeighboringCell = l ;
                            cells.at(i).neighbors.at(j).commonNeighbors.push_back( tmpAB ) ;
                            
                            CommonNeighbor tmpAC ;
                            tmpAC.cellIDCommonNeighbor = cellIDNeighbor ;
                            tmpAC.nghbrIDCell_CommonNeighbor = neighborIDCell ;
                            tmpAC.nghbrIDNeighboringCell_CommonNeighbor = l ;
                            tmpAC.nghbrIDCommonNeighbor_Cell = j ;
                            tmpAC.nghbrIDCommonNeighbor_NeighboringCell = nghbrIDNeighboringCell_Other ;
                            cells.at(i).neighbors.at(k).commonNeighbors.push_back( tmpAC ) ;
                            
                            CommonNeighbor tmpBA ;
                            tmpBA.cellIDCommonNeighbor = cellID_Other ;
                            tmpBA.nghbrIDCell_CommonNeighbor = nghbrIDNeighboringCell_Other ;
                            tmpBA.nghbrIDNeighboringCell_CommonNeighbor = nghbrIDCell_Other ;
                            tmpBA.nghbrIDCommonNeighbor_Cell = l ;
                            tmpBA.nghbrIDCommonNeighbor_NeighboringCell = k ;
                            cells.at(cellIDNeighbor).neighbors.at(neighborIDCell).commonNeighbors.push_back(tmpBA ) ;
                            
                            CommonNeighbor tmpBC ;
                            tmpBC.cellIDCommonNeighbor = i ;
                            tmpBC.nghbrIDCell_CommonNeighbor = j ;
                            tmpBC.nghbrIDNeighboringCell_CommonNeighbor = k ;
                            tmpBC.nghbrIDCommonNeighbor_Cell = neighborIDCell ;
                            tmpBC.nghbrIDCommonNeighbor_NeighboringCell = nghbrIDCell_Other ;
                            cells.at(cellIDNeighbor).neighbors.at(l).commonNeighbors.push_back(tmpBC) ;
                            
                            CommonNeighbor tmpCA ;
                            tmpCA.cellIDCommonNeighbor = cellIDNeighbor ;
                            tmpCA.nghbrIDCell_CommonNeighbor = l ;
                            tmpCA.nghbrIDNeighboringCell_CommonNeighbor = k ;
                            tmpCA.nghbrIDCommonNeighbor_NeighboringCell = j ;
                            tmpCA.nghbrIDCommonNeighbor_Cell = nghbrIDNeighboringCell_Other ;
                            
                            cells.at(cellID_Other).neighbors.at(nghbrIDCell_Other).commonNeighbors.push_back( tmpCA)  ;
                            
                            CommonNeighbor tmpCB ;
                            tmpCB.cellIDCommonNeighbor = i ;
                            tmpCB.nghbrIDCell_CommonNeighbor = k ;
                            tmpCB.nghbrIDNeighboringCell_CommonNeighbor = j ;
                            tmpCB.nghbrIDCommonNeighbor_Cell = nghbrIDCell_Other ;
                            tmpCB.nghbrIDCommonNeighbor_NeighboringCell = neighborIDCell ;
                            cells.at(cellID_Other).neighbors.at( nghbrIDNeighboringCell_Other).commonNeighbors.push_back( tmpCB ) ;
                            
                            
                        }
                    }
                }
            }
        }
    }

}
//---------------------------------------------------------------------------------------------
void SignalTissue::Print_CommonNeighbors()
{
    for (unsigned int i=0; i < cells.size(); i++)
    {
        for (unsigned int j=0; j< cells.at(i).neighbors.size(); j++)
        {
            cout<<i<<'\t'<<cells.at(i).neighbors.at(j).CellID_Neighbor<<" : "<<'\t' ;
            for (unsigned int k=0; k < cells.at(i).neighbors.at(j).commonNeighbors.size() ; k++)
            {
                cout<< cells.at(i).neighbors.at(j).commonNeighbors.at(k).cellIDCommonNeighbor <<'\t' ;
            }
            cout<<endl ;
        }
    }
}
//---------------------------------------------------------------------------------------------

void SignalTissue::Cal_Intersect()
{
    for (unsigned int i=0; i < cells.size(); i++)
    {
        for (unsigned int j=0; j< cells.at(i).neighbors.size(); j++)
        {
            int cellID_B = cells.at(i).neighbors.at(j).CellID_Neighbor ;
            int nghbrIDA_B = cells.at(i).neighbors.at(j).NeighborID_Cell ;
            if (i < cellID_B)
            {
                for (unsigned int k=0; k < cells.at(i).neighbors.at(j).commonNeighbors.size() ; k++)
                {
                    int cellID_C = cells.at(i).neighbors.at(j).commonNeighbors.at(k).cellIDCommonNeighbor ;
                    if (cellID_B < cellID_C )
                    {
                        int nghbrIDA_C = cells.at(i).neighbors.at(j).commonNeighbors.at(k).nghbrIDCell_CommonNeighbor ;
                        int nghbrIDC_A = cells.at(i).neighbors.at(j).commonNeighbors.at(k).nghbrIDCommonNeighbor_Cell ;
                        int nghbrIDB_C = cells.at(i).neighbors.at(j).commonNeighbors.at(k).nghbrIDNeighboringCell_CommonNeighbor ;
                        int nghbrIDC_B = cells.at(i).neighbors.at(j).commonNeighbors.at(k).nghbrIDCommonNeighbor_NeighboringCell ;
                        vector<int> indices ;
                        
                        vector<double> intfX_AB = cells.at(i).neighbors.at(j).intfX ;
                        vector<double> intfY_AB = cells.at(i).neighbors.at(j).intfY ;
                        
                        vector<double> intfX_AC = cells.at(i).neighbors.at(nghbrIDC_A).intfX ;
                        vector<double> intfY_AC = cells.at(i).neighbors.at(nghbrIDC_A).intfY ;
                        
                        vector<double> intfX_BC = cells.at(cellID_C).neighbors.at(nghbrIDB_C).intfX ;
                        vector<double> intfY_BC = cells.at(cellID_C).neighbors.at(nghbrIDB_C).intfY ;
                        
                        vector<vector<double> > AB_AC = Dist_VecToVec2D(intfX_AB, intfY_AB, intfX_AC, intfY_AC) ;
                        indices = Indices_MinMatrix(AB_AC) ;
                        
                        vector<vector<double> > tmp ;
                        tmp.push_back( Dist_PointToVec2D(intfX_AB.at(indices.at(0)), intfY_AB.at(indices.at(0)), intfX_BC, intfY_BC) );
                        indices.push_back(Indices_MinMatrix(tmp).at(1)) ;
                        
                        
                        double intrsctX = ( intfX_AB.at(indices.at(0))+ intfX_AC.at(indices.at(1))+ intfX_BC.at(indices.at(2)))/3.0 ;
                        double intrsctY = ( intfY_AB.at(indices.at(0))+ intfY_AC.at(indices.at(1))+ intfY_BC.at(indices.at(2)))/3.0 ;
                        
                        vector<double> tmpx {intfX_AB.at(indices.at(0)),intfX_AC.at(indices.at(1)),intfX_BC.at(indices.at(2))} ;
                        vector<double> tmpy {intfY_AB.at(indices.at(0)),intfY_AC.at(indices.at(1)),intfY_BC.at(indices.at(2))} ;
                        
                        vector<vector<double> > tmpDistVecToVec = Dist_VecToVec2D(tmpx, tmpy, tmpx, tmpy) ;
                        vector<double> tmpDist = *max_element(tmpDistVecToVec.begin(), tmpDistVecToVec.end()) ;
                        double maxDist = *max_element(tmpDist.begin(), tmpDist.end()) ;
                        
                        if ( maxDist < thres_intersect)
                        {
                            cells.at(i).neighbors.at(j).intersectX.push_back(intrsctX) ;
                            cells.at(i).neighbors.at(j).intersectY.push_back(intrsctY) ;
                            
                            cells.at(i).neighbors.at(nghbrIDC_A).intersectX.push_back(intrsctX) ;
                            cells.at(i).neighbors.at(nghbrIDC_A).intersectY.push_back(intrsctY) ;
                            
                            cells.at(cellID_B).neighbors.at(nghbrIDA_B).intersectX.push_back(intrsctX) ;
                            cells.at(cellID_B).neighbors.at(nghbrIDA_B).intersectY.push_back(intrsctY) ;
                            
                            cells.at(cellID_B).neighbors.at(nghbrIDC_B).intersectX.push_back(intrsctX) ;
                            cells.at(cellID_B).neighbors.at(nghbrIDC_B).intersectY.push_back(intrsctY) ;
                            
                            cells.at(cellID_C).neighbors.at(nghbrIDA_C).intersectX.push_back(intrsctX) ;
                            cells.at(cellID_C).neighbors.at(nghbrIDA_C).intersectY.push_back(intrsctY) ;
                            
                            cells.at(cellID_C).neighbors.at(nghbrIDB_C).intersectX.push_back(intrsctX) ;
                            cells.at(cellID_C).neighbors.at(nghbrIDB_C).intersectY.push_back(intrsctY) ;
                            
                            cells.at(i).neighbors.at(j).realIntersect.push_back( true ) ;
                            cells.at(i).neighbors.at(nghbrIDC_A).realIntersect.push_back(true) ;
                            cells.at(cellID_B).neighbors.at(nghbrIDA_B).realIntersect.push_back(true) ;
                            cells.at(cellID_B).neighbors.at(nghbrIDC_B).realIntersect.push_back(true) ;
                            cells.at(cellID_C).neighbors.at(nghbrIDA_C).realIntersect.push_back(true) ;
                            cells.at(cellID_C).neighbors.at(nghbrIDB_C).realIntersect.push_back(true) ;
                            
                        }
                        else
                        {
                            cells.at(i).neighbors.at(j).realIntersect.push_back( false) ;
                            cells.at(i).neighbors.at(nghbrIDC_A).realIntersect.push_back(false) ;
                            cells.at(cellID_B).neighbors.at(nghbrIDA_B).realIntersect.push_back(false) ;
                            cells.at(cellID_B).neighbors.at(nghbrIDC_B).realIntersect.push_back(false) ;
                            cells.at(cellID_C).neighbors.at(nghbrIDA_C).realIntersect.push_back(false) ;
                            cells.at(cellID_C).neighbors.at(nghbrIDB_C).realIntersect.push_back(false) ;
                            
                        }
                    }
                }
            }
        }
    }
}
//---------------------------------------------------------------------------------------------
void SignalTissue::Find_boundaries()
{
    vector<double> allNodesX ;
    vector<double> allNodesY ;
    for (unsigned int i = 0; i< cells.size(); i++)
    {
        for (unsigned int j =0; j < cells.at(i).neighbors.size(); j++)
        {
            int cellID_nghbr = cells.at(i).neighbors.at(j).CellID_Neighbor ;
            if (cells.at(i).neighbors.at(j).intersectX.size() < 2 || cells.at(i).noNeighboringNodesX.size() > 6 )
            {
                cells.at(i).boundary = true ;
                cells.at(cellID_nghbr).boundary = true ;
                if (cells.at(i).neighbors.at(j).intersectX.size() == 0 )
                {
                  vector<vector<double> > distTmp =  Dist_VecToVec2D(cells.at(i).neighbors.at(j).intfX,cells.at(i).neighbors.at(j).intfY,
                                                                    cells.at(i).neighbors.at(j).intfX, cells.at(i).neighbors.at(j).intfY ) ;
                    
                   vector<int> indecies =  Indices_MaxMatrix(distTmp) ;
                    
            //        cells.at(i).verticesX.push_back(cells.at(i).neighbors.at(j).intfX.at(indecies.at(0))) ;
            //        cells.at(i).verticesY.push_back(cells.at(i).neighbors.at(j).intfY.at(indecies.at(0))) ;
                    
            //        cells.at(i).verticesX.push_back(cells.at(i).neighbors.at(j).intfX.at(indecies.at(1))) ;
            //        cells.at(i).verticesY.push_back(cells.at(i).neighbors.at(j).intfY.at(indecies.at(1))) ;
                    
                    cells.at(i).newVertX.push_back(cells.at(i).neighbors.at(j).intfX.at(indecies.at(0))) ;
                    cells.at(i).newVertY.push_back(cells.at(i).neighbors.at(j).intfY.at(indecies.at(0))) ;
                    cells.at(i).newVertUpdateStatus.push_back(false) ;
                    
                    cells.at(i).newVertX.push_back(cells.at(i).neighbors.at(j).intfX.at(indecies.at(1))) ;
                    cells.at(i).newVertY.push_back(cells.at(i).neighbors.at(j).intfY.at(indecies.at(1))) ;
                    cells.at(i).newVertUpdateStatus.push_back(false) ;
                    
                }
                else if (cells.at(i).neighbors.at(j).intersectX.size() == 1)
                {
                    vector<double> distTmp = Dist_PointToVec2D(cells.at(i).neighbors.at(j).intersectX.at(0),
                                                               cells.at(i).neighbors.at(j).intersectY.at(0),
                                                               cells.at(i).neighbors.at(j).intfX,
                                                               cells.at(i).neighbors.at(j).intfY ) ;
                    
                    int indextmp = max_element(distTmp.begin(), distTmp.end() ) - distTmp.begin() ;
                    
            //        cells.at(i).verticesX.push_back( cells.at(i).neighbors.at(j).intfX.at(indextmp)) ;
            //        cells.at(i).verticesY.push_back( cells.at(i).neighbors.at(j).intfY.at(indextmp)) ;
                    
                    cells.at(i).newVertX.push_back( cells.at(i).neighbors.at(j).intfX.at(indextmp)) ;
                    cells.at(i).newVertY.push_back(cells.at(i).neighbors.at(j).intfY.at(indextmp)) ;
                    cells.at(i).newVertUpdateStatus.push_back(false) ;
                }
            }
        }
        
    }
}
//---------------------------------------------------------------------------------------------
void SignalTissue::Count_IntersectPoints()
{
    for (unsigned int i = 0; i< cells.size(); i++)
    {
        for (unsigned int j = 0; j < cells.at(i).neighbors.size(); j++)
        {
           cout<< cells.at(i).neighbors.at(j).intersectX.size()<< endl ;
        }
    }
    
}
//---------------------------------------------------------------------------------------------
void SignalTissue::Cal_AllCellVertices()
{
    for (unsigned int i =0; i < cells.size(); i++)
    {
        cells.at(i).Cal_Vertices() ;
    }
}
//---------------------------------------------------------------------------------------------
void SignalTissue::ParaViewVertices ()
{
    if (writeVtk == false)
    {
        return ;
    }
    vector<double> allNodesX ;
    vector<double> allNodesY ;
    vector<int> links ;         // link to the other vertices
    vector<int > links2 ;       // link to the centroid of the cell
    int verticesSize ;
    for (unsigned int i=0 ; i< cells.size() ; i++)
    {
        for (unsigned int j=0 ; j< cells.at(i).verticesX.size() ; j++)
        {
            allNodesX.push_back( cells.at(i).verticesX.at(j))   ;
            allNodesY.push_back( cells.at(i).verticesY.at(j))   ;
            links.push_back(cells.at(i).connections.at(j) ) ;
        }
    }
    verticesSize = allNodesX.size() ;       //centroids start from here in the vector
    for (unsigned int i=0 ; i< cells.size() ; i++)
    {
        allNodesX.push_back(cells.at(i).centroid.at(0) )   ;
        allNodesY.push_back(cells.at(i).centroid.at(1) )   ;
        for (unsigned int j=0 ; j< cells.at(i).verticesX.size() ; j++)
        {
            // (verticesSize + i) is equal to index of centroid of the i_th cell
            links2.push_back(verticesSize + i ) ;
        }
    }
    string vtkFileName = "Vertices"+ to_string(frameIndex)+ ".vtk" ;
    ofstream VerticesOut;
    VerticesOut.open(vtkFileName.c_str());
    VerticesOut<< "# vtk DataFile Version 3.0" << endl;
    VerticesOut<< "Result for paraview 2d code" << endl;
    VerticesOut << "ASCII" << endl;
    VerticesOut << "DATASET UNSTRUCTURED_GRID" << endl;
    VerticesOut << "POINTS " << allNodesX.size()   << " float" << endl;
    
    for (unsigned int i = 0; i < allNodesX.size(); i++)
    {
        
        VerticesOut << allNodesX.at(i) << " " << allNodesY.at(i) << " " << 0.0 << endl;
    }
    VerticesOut<< endl;
    
    VerticesOut<< "CELLS " << links.size()+ links2.size() << " " << 3 *( links.size()+ links2.size() )<< endl;
    
    for (unsigned int i = 0; i < (links.size()); i++)           //number of connections per node
    {
        VerticesOut << 2 << " " << i << " " << links.at(i) << endl;
        VerticesOut << 2 << " " << i << " " << links2.at(i) << endl;
        
    }
    
    VerticesOut << "CELL_TYPES " << links.size()+ links2.size()<< endl;             //connection type
    for (unsigned int i = 0; i < links.size()+ links2.size(); i++) {
        VerticesOut << "3" << endl;
    }
    
    VerticesOut << "POINT_DATA "<<allNodesX.size() <<endl ;
    VerticesOut << "SCALARS Cell_ID " << "float"<< endl;
    VerticesOut << "LOOKUP_TABLE " << "default"<< endl;
    for (unsigned int i = 0; i < cells.size() ; i++)
    {
        for ( uint j=0; j < cells.at(i).verticesX.size() ; j++)
        {
            VerticesOut<< i <<endl ;
        }
    }
    for (uint i = 0; i < cells.size() ; i++)
    {
        VerticesOut<< i <<endl ;
    }
    
    VerticesOut.close();
}
//---------------------------------------------------------------------------------------------

void SignalTissue::ParaViewTissue ()
{
    if (writeVtk == false)
    {
        return ;
    }
    vector<double> allNodesX ;
    vector<double> allNodesY ;
    for (unsigned int i=0 ; i< cells.size() ; i++)
    {
        for (unsigned int j=0 ; j< cells.at(i).nodesXNew.size() ; j++)
        {
            allNodesX.push_back(cells.at(i).nodesXNew.at(j))   ;
            allNodesY.push_back(cells.at(i).nodesYNew.at(j))   ;
        }
    }
    int size =static_cast<int> (allNodesX.size() );
    string vtkFileName = "Tissue"+ to_string(frameIndex)+ ".vtk" ;
    ofstream TissueOut;
    TissueOut.open(vtkFileName.c_str());
    TissueOut<< "# vtk DataFile Version 3.0" << endl;
    TissueOut<< "Result for paraview 2d code" << endl;
    TissueOut << "ASCII" << endl;
    TissueOut << "DATASET UNSTRUCTURED_GRID" << endl;
    TissueOut << "POINTS " << size   << " float" << endl;
    
    for (uint i = 0; i < size ; i++)
    {
        
        TissueOut << allNodesX.at(i) << " " << allNodesY.at(i) << " " << 0.0 << endl;
    }
    TissueOut<< endl;
    TissueOut<< "CELLS " << size << " " << 3 *(size )<< endl;
    
    for (uint i = 0; i < size ; i++)           //number of connections per node
    {
        
        TissueOut << 2 << " " << i << " " << (i+1) % size << endl;
        
    }
    
    TissueOut << "CELL_TYPES " << size << endl;             //connection type
    for (uint i = 0; i < size ; i++) {
        TissueOut << "3" << endl;
    }
    
    TissueOut << "POINT_DATA "<< size <<endl ;
    TissueOut << "SCALARS Polygon " << "float"<< endl;
    TissueOut << "LOOKUP_TABLE " << "default"<< endl;
    for (uint i = 0; i < cells.size() ; i++)
    {
        for ( uint j=0; j< cells.at(i).nodesXNew.size()  ; j++)
        {
            TissueOut<< cells.at(i).verticesX.size() <<endl ;
         //   TissueOut<< cells.at(i).layer<<endl ;
        }
    }
    
    TissueOut.close();
}
//---------------------------------------------------------------------------------------------
void SignalTissue::Cal_AllCellConnections()
{
     int tmpShift = 0 ;
    int counter = 0 ;
    for (unsigned int i =0; i < cells.size(); i++)
    {
        tmpShift = counter ;
        for (unsigned int j=0 ; j< cells.at(i).verticesX.size(); j++)
        {
            cells.at(i).connections.push_back( (j+1) % cells.at(i).verticesX.size() + tmpShift ) ;
            counter += 1 ;
        }
    }
}
//---------------------------------------------------------------------------------------------
void SignalTissue:: Test()
{
    for (unsigned int i =0; i < cells.size(); i++)
    {
        cout<<cells.at(i).verticesX.size()<<endl ;
    }
}
//---------------------------------------------------------------------------------------------
void SignalTissue:: Add_NewVerticesToBoundaryEdges()
{
    for (unsigned int i =0; i < cells.size(); i++)
    {
        cells.at(i).Add_BoundaryVertice3 () ;
    }
}

//---------------------------------------------------------------------------------------------
void SignalTissue:: Refine_VerticesInBoundaryCells ()
{
    // Add vertices for nodes that are shared only between two cells but not three.
    for (unsigned int i =0; i < cells.size(); i++)
    {
        if (cells.at(i).boundary )
        {
            for (unsigned int j =0; j < cells.at(i).newVertX.size(); j++)
            {
                if (cells.at(i).newVertUpdateStatus.at(j) == false)
                {
                    double x1 = cells.at(i).newVertX.at(j) ;
                    double y1 = cells.at(i).newVertY.at(j) ;
                    vector<double> pointsX ;
                    vector<double> pointsY ;
                    vector<int > cellIDs ;
                    vector<int> newvertID ;
                    for (int k =0 ; k < cells.size() ; k++)
                    {
                        if (cells.at(k).boundary)
                        {
                            for (unsigned int l = 0; l < cells.at(k).newVertX.size(); l++)
                            {
                                double tmpDist = Dist2D(x1, y1, cells.at(k).newVertX.at(l), cells.at(k).newVertY.at(l)) ;
                                if (tmpDist < thres_corners )
                                {
                                    pointsX.push_back(cells.at(k).newVertX.at(l)) ;
                                    pointsY.push_back(cells.at(k).newVertY.at(l)) ;
                                    cellIDs.push_back(k) ;
                                    newvertID.push_back(l) ;
                                }
                                
                            }
                        }
                        
                    }
                    if (pointsX.size() > 4 || pointsX.size()==2)        // temporary works
                    {
                        double newvertX =  accumulate(pointsX.begin(), pointsX.end(), 0.0)/pointsX.size() ;
                        double newvertY =  accumulate(pointsY.begin(), pointsY.end(), 0.0)/pointsY.size() ;
                        
                        for (unsigned int n=0; n < cellIDs.size(); n++)
                        {
                            cells.at(cellIDs.at(n)).newVertUpdateStatus.at(newvertID.at(n)) = true ;
                        }
                        sort(cellIDs.begin(), cellIDs.end()) ;
                        cellIDs.erase( unique(cellIDs.begin(),cellIDs.end()) , cellIDs.end() ) ;
                        for (unsigned int n=0 ; n< cellIDs.size(); n++)
                        {
                            cells.at(cellIDs.at(n)).verticesX.push_back(newvertX) ;
                            cells.at(cellIDs.at(n)).verticesY.push_back(newvertY) ;
                            
                        }
                    }
                    else
                    {
                        cells.at(i).verticesX.push_back(x1) ;
                        cells.at(i).verticesY.push_back(y1) ;
                    }
                    
                }
            }
            
        }
    }
}
//---------------------------------------------------------------------------------------------

void SignalTissue::Cyclic4Correction()
{
    for (unsigned int i =0 ; i< cells.size(); i++)
    {
     //   int verticesSize = cells.at(i).verticesX.size() ;
        for (unsigned int j =0; j < cells.at(i).verticesX.size() ; j++)
        {
            vector<double> pointsX ;
            vector<double> pointsY ;
            vector<int > cellIDs ;
            vector<int> vertID ;
            
            for (unsigned int k =0; k < cells.at(i).cyclic4.size(); k++)
            {
                int cellID = cells.at(i).cyclic4.at(k) ;
              //  if (i < cellID)
                {
                    for (unsigned int l=0; l < cells.at(cellID).verticesX.size(); l++)
                    {
                        double x1 = cells.at(i).verticesX.at(j) ;
                        double y1 = cells.at(i).verticesY.at(j) ;
                        
                        double tmpDist = Dist2D(x1, y1, cells.at(cellID).verticesX.at(l), cells.at(cellID).verticesY.at(l)) ;
                        if (tmpDist < thres_cyclic4 )
                        {
                            pointsX.push_back(cells.at(cellID).verticesX.at(l)) ;
                            pointsY.push_back(cells.at(cellID).verticesY.at(l)) ;
                            cellIDs.push_back(cellID) ;
                            vertID.push_back(l) ;
                        }
                    }
                }
            }
            if (pointsX.size() > 0)
            {
                double newvertX =  accumulate(pointsX.begin(), pointsX.end(), 0.0)/pointsX.size() ;
                double newvertY =  accumulate(pointsY.begin(), pointsY.end(), 0.0)/pointsY.size() ;
                for (long n = vertID.size()-1 ; n >= 0 ; n--)
                {
                    cells.at(cellIDs.at(n)).verticesX.erase( cells.at(cellIDs.at(n)).verticesX.begin() + vertID.at(n) ) ;
                    cells.at(cellIDs.at(n)).verticesY.erase( cells.at(cellIDs.at(n)).verticesY.begin() + vertID.at(n) ) ;
                }
                
                sort(cellIDs.begin(), cellIDs.end()) ;
                cellIDs.erase( unique(cellIDs.begin(),cellIDs.end()) , cellIDs.end() ) ;
                for (unsigned int n=0 ; n< cellIDs.size(); n++)
                {
                    cells.at(cellIDs.at(n)).verticesX.push_back(newvertX) ;
                    cells.at(cellIDs.at(n)).verticesY.push_back(newvertY) ;
                    
                }
            }
                                
        }
    }
    
}

//---------------------------------------------------------------------------------------------
void SignalTissue:: Find_Cyclic4()
{
    for (unsigned int i =0; i < cells.size(); i++)
    {
        for (unsigned int j =0; j< cells.at(i).neighbors.size(); j++)
        {
            
            int cellID_B = cells.at(i).neighbors.at(j).CellID_Neighbor ;
            if ( i < cellID_B)
            {
                for (unsigned int k=0; k< cells.at(i).neighbors.at(j).commonNeighbors.size(); k++)
                {
                    for (unsigned int l = k+1; l < cells.at(i).neighbors.at(j).commonNeighbors.size() ; l++)
                    {
                        int cellID_C = cells.at(i).neighbors.at(j).commonNeighbors.at(k).cellIDCommonNeighbor ;
                        int cellID_D = cells.at(i).neighbors.at(j).commonNeighbors.at(l).cellIDCommonNeighbor ;
                        if (cellID_B < cellID_C && cellID_C < cellID_D )
                            {
                            
                            for (unsigned int m =0; m< cells.at(cellID_C).neighbors.size(); m++)
                            {
                                if (cells.at(cellID_C).neighbors.at(m).CellID_Neighbor == cellID_D )
                                {
                                    cells.at(i).cyclic4.push_back(i) ;
                                    cells.at(i).cyclic4.push_back(cellID_B) ;
                                    cells.at(i).cyclic4.push_back(cellID_C) ;
                                    cells.at(i).cyclic4.push_back(cellID_D) ;
                                    
                                    cells.at(cellID_B).cyclic4.push_back(i) ;
                                    cells.at(cellID_B).cyclic4.push_back(cellID_B) ;
                                    cells.at(cellID_B).cyclic4.push_back(cellID_C) ;
                                    cells.at(cellID_B).cyclic4.push_back(cellID_D) ;
                                    
                                    cells.at(cellID_C).cyclic4.push_back(i) ;
                                    cells.at(cellID_C).cyclic4.push_back(cellID_B) ;
                                    cells.at(cellID_C).cyclic4.push_back(cellID_C) ;
                                    cells.at(cellID_C).cyclic4.push_back(cellID_D) ;
                                    
                                    cells.at(cellID_D).cyclic4.push_back(i) ;
                                    cells.at(cellID_D).cyclic4.push_back(cellID_B) ;
                                    cells.at(cellID_D).cyclic4.push_back(cellID_C) ;
                                    cells.at(cellID_D).cyclic4.push_back(cellID_D) ;
                                    
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for (unsigned int i =0; i < cells.size(); i++)
    {
        sort(cells.at(i).cyclic4.begin(), cells.at(i).cyclic4.end()) ;
        cells.at(i).cyclic4.erase( unique(cells.at(i).cyclic4.begin(),cells.at(i).cyclic4.end()) , cells.at(i).cyclic4.end() ) ;
    }
}

//---------------------------------------------------------------------------------------------

void SignalTissue::SortVertices()
{
    for (unsigned int i =0; i < cells.size(); i++)
    {
        cells.at(i).SortCCW(); // verticesX ,verticesY ) ;
    }
}
//---------------------------------------------------------------------------------------------
void SignalTissue::AllCell_RefineNoBoundary ()
{
    for (unsigned int i=0; i< cells.size(); i++)
    {
        cells.at(i).Refine_NoBoundary() ;
        //  cells.at(i).Refine_NodeXNew() ; 
    }
}
//---------------------------------------------------------------------------------------------
void SignalTissue::ParaViewBoundary ()
{
    if (writeVtk == false)
    {
        return ;
    }
    vector<double> allNodesX ;
    vector<double> allNodesY ;
    for (unsigned int i=0 ; i< cells.size() ; i++)
    {
        for (unsigned int j=0 ; j< cells.at(i).noNeighboringNodesX.size() ; j++)
        {
            allNodesX.push_back(cells.at(i).noNeighboringNodesX.at(j))   ;
            allNodesY.push_back(cells.at(i).noNeighboringNodesY.at(j))   ;
        }
    }
    int size =static_cast<int> (allNodesX.size() );
    string vtkFileName = "Boundary"+ to_string(frameIndex)+ ".vtk" ;
    ofstream BoundaryOut;
    BoundaryOut.open(vtkFileName.c_str());
    BoundaryOut<< "# vtk DataFile Version 3.0" << endl;
    BoundaryOut<< "Result for paraview 2d code" << endl;
    BoundaryOut << "ASCII" << endl;
    BoundaryOut << "DATASET UNSTRUCTURED_GRID" << endl;
    BoundaryOut << "POINTS " << size   << " float" << endl;
    
    for (uint i = 0; i < size ; i++)
    {
        
        BoundaryOut << allNodesX.at(i) << " " << allNodesY.at(i) << " " << 0.0 << endl;
    }
    BoundaryOut<< endl;
    BoundaryOut.close();
}
//---------------------------------------------------------------------------------------------
void SignalTissue::ParaViewInitialConfiguration ()
{
    if (writeVtk == false)
    {
        return ;
    }
    vector<double> allNodesX ;
    vector<double> allNodesY ;
    for (unsigned int i=0 ; i< cells.size() ; i++)
    {
        for (unsigned int j=0 ; j < cells.at(i).nodesX.size() ; j++)
        {
            allNodesX.push_back(cells.at(i).nodesX.at(j))   ;
            allNodesY.push_back(cells.at(i).nodesY.at(j))   ;
        }
    }
    int size =static_cast<int> (allNodesX.size() );
    string vtkFileName = "Initial"+ to_string(frameIndex)+ ".vtk" ;
    ofstream InitialOut;
    InitialOut.open(vtkFileName.c_str());
    InitialOut<< "# vtk DataFile Version 3.0" << endl;
    InitialOut<< "Result for paraview 2d code" << endl;
    InitialOut << "ASCII" << endl;
    InitialOut << "DATASET UNSTRUCTURED_GRID" << endl;
    InitialOut << "POINTS " << size   << " float" << endl;
    
    for (unsigned int i = 0; i < size ; i++)
    {
        
        InitialOut << allNodesX.at(i) << " " << allNodesY.at(i) << " " << 0.0 << endl;
    }
    
    InitialOut<< endl;
    InitialOut.close();
}
//---------------------------------------------------------------------------------------------
void SignalTissue::Refine_CurvedInterface ()
{
    for (unsigned int i = 0; i< cells.size(); i++)
    {
        for (unsigned int j=0; j< cells.at(i).neighbors.size(); j++)
        {
            if (cells.at(i).neighbors.at(j).curvedInterface == false)
            {
                int cellID_Neighbor = cells.at(i).neighbors.at(j).CellID_Neighbor ;
                int neighborID_Cell = cells.at(i).neighbors.at(j).NeighborID_Cell ;
                vector<double> intfX = cells.at(i).neighbors.at(j).intfX ;
                vector<double> intfY = cells.at(i).neighbors.at(j).intfY ;
                pair<double, double> point ;
                vector<double> tmpX ;
                vector<double> tmpY ;
                //find the first point of the sequences and connect the whole curve to the nearest vertices
                if (cells.at(i).neighbors.at(j).intersectX.size() > 0)
                {
                    intfX.insert(intfX.end(), cells.at(i).neighbors.at(j).intersectX.begin(), cells.at(i).neighbors.at(j).intersectX.end() ) ;
                    intfY.insert(intfY.end(), cells.at(i).neighbors.at(j).intersectY.begin(), cells.at(i).neighbors.at(j).intersectY.end() ) ;
                    point.first = intfX.back() ;
                    point.second = intfY.back() ;
                }
                else
                {
                    vector<vector<double> > pairwiseDist =  Dist_VecToVec2D(intfX , intfY, intfX , intfY ) ;
                    vector<int> indecies =  Indices_MaxMatrix(pairwiseDist) ;
                    point.first = intfX.at( indecies.at(0) ) ;
                    point.second = intfY.at(indecies.at(0) ) ;
                }
                vector<double> tmpDist ;
                tmpDist = Dist_PointToVec2D(point.first, point.second, intfX, intfY) ;
                vector<int > idx(tmpDist.size());
                iota(idx.begin(), idx.end(), 0);
                sort(idx.begin(), idx.end(),
                     [&tmpDist](int i1, int i2) {return tmpDist[i1] < tmpDist[i2];});
                for (unsigned int k = 0; k < tmpDist.size(); k++)
                {
                    tmpX.push_back(intfX.at(idx.at(k) ) ) ;
                    tmpY.push_back(intfY.at(idx.at(k) ) ) ;
                }
                double area = 0.0 ;
                
                vector<double > cntx(tmpX.size(), cells.at(i).centroid.front() ) ;
                vector<double > cnty(tmpY.size(), cells.at(i).centroid.back() ) ;
                vector<double > distToCntrX ;
                distToCntrX.resize(tmpX.size() ) ;
                vector<double > distToCntrY ;
                distToCntrY.resize(tmpY.size() ) ;
                transform(tmpX.begin(), tmpX.end(), cntx.begin(), distToCntrX.begin(), linearConfig(1.0 ,-1.0) ) ;
                transform(tmpY.begin(), tmpY.end(), cnty.begin(), distToCntrY.begin(), linearConfig(1.0 ,-1.0) ) ;
                vector<double> tetta ;
                vector<double> tmpMag ;
                
                for (unsigned int k=0 ; k<distToCntrX.size(); k++)
                {
                    tetta.push_back(atan2(distToCntrY.at(k), distToCntrX.at(k) ) ) ;
                    tmpMag.push_back(MagnitudeVec(distToCntrX.at(k), distToCntrY.at(k))) ;
                }
                vector<double> dTetta ;
                for (unsigned int k=0 ; k<distToCntrX.size()-1; k++)
                {
                    dTetta.push_back( tetta.at(k+1) - tetta.at(k) ) ;
                    area += 0.5 * tmpMag.at(k+1) * tmpMag.at(k) * abs( sin( dTetta.back() ) ) ;
                }
                double refArea = 0.5 * tmpMag.front() * tmpMag.back() * abs( sin( tetta.back()-tetta.front() ) ) ;
                //   if (abs( (area - refArea )/ area ) > 0.3 )
                if ( (area - refArea )/ area > 0.5 )        // do it based on cells which are bigger than current vertices
                {
                    vector<double> a ;
                    vector<double> dA ;
                    for (unsigned int k=0; k< tetta.size(); k++)
                    {
                        a.push_back(0.5 * tmpMag.at(k) * ( abs( sin( tetta.at(k)- tetta.front() ) ) * tmpMag.front() +
                                                          abs( sin( tetta.at(k)- tetta.back()  ) ) * tmpMag.back() ) ) ;
                        dA.push_back(abs(a.back() - area ) ) ;
                    }
                    int idNode = min_element(dA.begin(), dA.end() )- dA.begin() ;
                    refArea = a.at(idNode) ;
                    
                    //       cells.at(i).neighbors.at(j).curvedInterface = true ;
                    cells.at(i).verticesX.push_back(tmpX.at(idNode) ) ;
                    cells.at(i).verticesY.push_back(tmpY.at(idNode) ) ;
                    
                    //       cells.at(cellID_Neighbor).neighbors.at(neighborID_Cell ).curvedInterface = true ;
                    cells.at(cellID_Neighbor).verticesX.push_back(tmpX.at(idNode)) ;
                    cells.at(cellID_Neighbor).verticesY.push_back(tmpY.at(idNode)) ;
                    
                }
            }
        }
    }
}

//---------------------------------------------------------------------------------------------
void SignalTissue::Print_VeritcesSize()
{
    for (unsigned int i=0; i<cells.size(); i++)
    {
        cout<<"Cell "<<i<<" :"<<'\t'<<cells.at(i).verticesX.size()<<endl ;
    }
}





//--------------------------------------Meshes----------------=--------------------------------
//---------------------------------------------------------------------------------------------
void SignalTissue::Find_AllMeshes()
{
    for (unsigned int i =0 ; i < cells.size(); i++)
    {
        cells.at(i).Find_Mesh() ;
    }
}
//---------------------------------------------------------------------------------------------
void SignalTissue::Cal_AreaOfTissue()
{
    areaTissue = 0.0 ;
    for (unsigned int i =0; i<cells.size(); i++)
    {
        cells.at(i).areaCell = 0.0 ;
        for (unsigned int j=0; j<cells.at(i).meshes.size(); j++)
        {
           cells.at(i).areaCell += cells.at(i).meshes.at(j).Cal_MeshArea() ;
            
        }
        areaTissue += cells.at(i).areaCell ;
    }
    TissueRadius = sqrt(areaTissue/ pi ) ;
}
//---------------------------------------------------------------------------------------------
void SignalTissue::Cal_AllSelfDiffusion()
{
    for (unsigned int i =0 ; i < cells.size(); i++)
    {
        cells.at(i).Self_Diffusion() ;
    }
}
//---------------------------------------------------------------------------------------------
void SignalTissue::Find_IntercellularMeshConnection()
{
    for (unsigned int i =0; i < cells.size(); i++)
    {
        for (unsigned int k =0 ; k< cells.at(i).meshes.size(); k++)
        {
            int counter = 0 ;
            for (unsigned int j =0; j < cells.at(i).neighbors.size(); j++)
            {
                int cell2 = cells.at(i).neighbors.at(j).CellID_Neighbor ;
                for (unsigned int l = 0; l < cells.at(cell2).meshes.size() ; l++)
                {
                    double smallValue = 0.0001 ;
                    counter = 0 ;
                    
                    if ( Dist2D(cells.at(i).meshes.at(k).triangleX.at(1), cells.at(i).meshes.at(k).triangleY.at(1),
                                cells.at(cell2).meshes.at(l).triangleX.at(1), cells.at(cell2).meshes.at(l).triangleY.at(1)) < smallValue
                        || Dist2D(cells.at(i).meshes.at(k).triangleX.at(1), cells.at(i).meshes.at(k).triangleY.at(1),
                                  cells.at(cell2).meshes.at(l).triangleX.at(2), cells.at(cell2).meshes.at(l).triangleY.at(2)) < smallValue )
                    {
                        counter +=1 ;
                    }
                    if ( Dist2D(cells.at(i).meshes.at(k).triangleX.at(2), cells.at(i).meshes.at(k).triangleY.at(2),
                                cells.at(cell2).meshes.at(l).triangleX.at(1), cells.at(cell2).meshes.at(l).triangleY.at(1)) < smallValue
                        || Dist2D(cells.at(i).meshes.at(k).triangleX.at(2), cells.at(i).meshes.at(k).triangleY.at(2),
                                  cells.at(cell2).meshes.at(l).triangleX.at(2), cells.at(cell2).meshes.at(l).triangleY.at(2)) < smallValue )
                    {
                        counter +=1 ;
                    }
                    if (counter ==2)
                    {
                        cells.at(i).meshes.at(k).connection.first = cell2 ;
                        cells.at(i).meshes.at(k).connection.second = l ;
                        double length =Dist2D(cells.at(i).meshes.at(k).center.at(0),cells.at(i).meshes.at(k).center.at(1),
                                              cells.at(cell2).meshes.at(l).center.at(0), cells.at(cell2).meshes.at(l).center.at(1)) ;
                        
                        cells.at(i).meshes.at(k).length.at(2) = length ;
                        
                        cells.at(cell2).meshes.at(l).connection.first = i ;
                        cells.at(cell2).meshes.at(l).connection.second = k ;
                        cells.at(cell2).meshes.at(l).length.at(2) = length ;
                        
                        break ;
                    }
                }
                if (counter ==2)    break ;
                
            }
        }
    }
}

//---------------------------------------------------------------------------------------------

void SignalTissue::IntercellularDiffusion()
{
    double dInter = 1.0 ;
    for (unsigned int i =0 ; i < cells.size(); i++)
    {
        for (unsigned int j =0; j < cells.at(i).meshes.size(); j++)
        {
            int cellID = cells.at(i).meshes.at(j).connection.first ;
            int meshID = cells.at(i).meshes.at(j).connection.second ;
            //cellID== -1 means the mesh is a boundary mesh (internal connections only )
            if (cellID == -1)
            {
                continue ;
            }
            double tmpChannel = cells.at(i).meshes.at(j).channel.at(2) ;
            double length = cells.at(i).meshes.at(j).length.at(2) ;
            double deltaU = cells.at(cellID).meshes.at(meshID).u1 - cells.at(i).meshes.at(j).u1 ;
            // meshSize squared
            cells.at(i).meshes.at(j).flux.push_back(dInter * tmpChannel * deltaU / (length * length) ) ;
            
        }
    }
}
//---------------------------------------------------------------------------------------------
void SignalTissue::EulerMethod()
{
    bool state = false ;
    int l = 0 ;
    while (state==false)
    {
        Cal_AllSelfDiffusion() ;
        IntercellularDiffusion() ;
      //  if (l%1000==0) cout<<l/1000<<endl ;
      //  if (l%1000==0) ParaViewMesh(l/1000) ;
        for (unsigned int i = 0; i < cells.size(); i++)
        {
            for (unsigned int j =0; j < cells.at(i).meshes.size(); j++)
            {
                cells.at(i).meshes.at(j).Euler() ;
            }
        }
        
        state = true ;
        for (unsigned int i = 0; i < cells.size(); i++)
        {
            for (unsigned int j =0; j < cells.at(i).meshes.size(); j++)
            {
                 if( abs(cells.at(i).meshes.at(j).u2- cells.at(i).meshes.at(j).u1)/(cells.at(i).meshes.at(j).u1 + 0.000001) > 0.000001)
                 {
                     state = false ;
                     break ;
                 }
                 
            }
            if (state == false) break ;
        }
        for (unsigned int i = 0; i < cells.size(); i++)
        {
            for (unsigned int j =0; j < cells.at(i).meshes.size(); j++)
            {
                cells.at(i).meshes.at(j).UpdateU() ;
            }
        }
        
        l++ ;
    }
     ParaViewMesh(frameIndex) ;
    /*
    double value = 0 ;
    for (unsigned int j =0 ; j < cells.size();j++)
    {
        for (unsigned int i =0; i< cells.at(j).meshes.size(); i++)
        {
            value += cells.at(j).meshes.at(i).u1 ;
        }
    }

    cout<< "value is "<<value<<endl ;
     */
    cout<< "number of steps needed is " << l <<endl ;
}
//---------------------------------------------------------------------------------------------
void SignalTissue::EulerMethod2 ()
{
    
    for (int l =0 ; l <= 100000; l++)
    {
        Cal_AllSelfDiffusion() ;
        IntercellularDiffusion() ;
        if (l%1000==0) cout<<l/1000<<endl ;
        if (l%1000==0)
        {
            ParaViewMesh(l) ;
        }
        
        for (unsigned int i = 0; i < cells.size(); i++)
        {
            for (unsigned int j =0; j < cells.at(i).meshes.size(); j++)
            {
                cells.at(i).meshes.at(j).Euler() ;
           //     cells.at(i).meshes.at(j).FullModel_Euler() ;
                cells.at(i).meshes.at(j).UpdateU() ;
            }
        }
    }
}
//---------------------------------------------------------------------------------------------
void SignalTissue::Find_SecretingCell()
{
    double cntX = tissueCenter.at(0) ;
    double cntY = tissueCenter.at(1) ;
    
    double minDis = Dist2D(cntX, cntY , cells.at(0).centroid.at(0), cells.at(0).centroid.at(1)) ;
    int cellID = 0 ;
    double dis ;
    for (unsigned int i = 1; i< cells.size(); i++)
    {
        dis = Dist2D(cntX, cntY, cells.at(i).centroid.at(0), cells.at(i).centroid.at(1)) ;
        if (dis < minDis)
        {
            minDis = dis ;
            cellID = i ;
        }
    }
    cout <<"cellID source " << cellID <<endl ;
    for (unsigned int j = 0; j < cells.at(cellID).meshes.size(); j++)
    {
        cells.at(cellID).meshes.at(j).c = 1.0 ;
    }
}
//---------------------------------------------------------------------------------------------
void SignalTissue::ParaViewMesh(int number)
{
    if (writeVtk == false)
    {
        return ;
    }
    vector<double> allNodesX ;
    vector<double> allNodesY ;
    int numberOfCells = 0 ;
    for (unsigned int i=0 ; i< cells.size() ; i++)
    {
        for (unsigned int j=0 ; j< cells.at(i).meshes.size() ; j++)
        {
            numberOfCells += 1 ;
            for (unsigned int k =0; k< cells.at(i).meshes.at(j).triangleX.size(); k++)
            {
                allNodesX.push_back( cells.at(i).meshes.at(j).triangleX.at(k))   ;
                allNodesY.push_back( cells.at(i).meshes.at(j).triangleY.at(k))   ;
            }
        }
    }
    
    string vtkFileName = "Mesh"+ to_string(number)+ ".vtk" ;
    ofstream MeshOut;
    MeshOut.open(vtkFileName.c_str());
    MeshOut<< "# vtk DataFile Version 3.0" << endl;
    MeshOut<< "Result for paraview 2d code" << endl;
    MeshOut << "ASCII" << endl;
    MeshOut << "DATASET UNSTRUCTURED_GRID" << endl;
    MeshOut << "POINTS " << allNodesX.size()   << " float" << endl;
    
    for (uint i = 0; i < allNodesX.size(); i++)
    {
        
        MeshOut << allNodesX.at(i) << " " << allNodesY.at(i) << " " << 0.0 << endl;
    }
    MeshOut<< endl;
    MeshOut<< "CELLS " << numberOfCells << " " << 4 *( numberOfCells )<< endl;
    
    for (int i=0 ; i< numberOfCells ; i++)
    {
            MeshOut << 3 << " " << 3*i << " " << 3*i +1  << " " << 3*i +2 <<endl;
    }
    
    MeshOut << "CELL_TYPES " << numberOfCells << endl;             //connection type
    for (uint i = 0; i < numberOfCells ; i++) {
        MeshOut << "5" << endl;
    }
    if (equationsType == simpleODE)
    {
        MeshOut << "POINT_DATA "<<allNodesX.size() <<endl ;
        MeshOut << "SCALARS WUS " << "float"<< endl;
        MeshOut << "LOOKUP_TABLE " << "default"<< endl;
        for (uint i = 0; i < cells.size() ; i++)
        {
            for ( uint j=0; j < cells.at(i).meshes.size() ; j++)
            {
                MeshOut << cells.at(i).meshes.at(j).u1 << " " <<cells.at(i).meshes.at(j).u1 << " "<<cells.at(i).meshes.at(j).u1 <<endl ;
            }
        }
    }
    if (equationsType == fullModel)
    {
        MeshOut << "POINT_DATA "<<allNodesX.size() <<endl ;
        MeshOut << "SCALARS WUS_mRNA " << "float"<< endl;
        MeshOut << "LOOKUP_TABLE " << "default"<< endl;
        for (uint i = 0; i < cells.size() ; i++)
        {
            for ( uint j=0; j < cells.at(i).meshes.size() ; j++)
            {
                   MeshOut << cells.at(i).meshes.at(j).concentrations2.at(0)<<" " <<cells.at(i).meshes.at(j).concentrations2.at(0) << " "<<cells.at(i).meshes.at(j).concentrations2.at(0) <<endl ;
            }
        }
        
        MeshOut << "SCALARS WUS_nucleus " << "float"<< endl;
        MeshOut << "LOOKUP_TABLE " << "default"<< endl;
        for (uint i = 0; i < cells.size() ; i++)
        {
            for ( uint j=0; j < cells.at(i).meshes.size() ; j++)
            {
                MeshOut << cells.at(i).meshes.at(j).concentrations2.at(1)<<" " <<cells.at(i).meshes.at(j).concentrations2.at(1) << " "<<cells.at(i).meshes.at(j).concentrations2.at(1) <<endl ;
            }
        }

        MeshOut << "SCALARS WUS_cytoplasm " << "float"<< endl;
        MeshOut << "LOOKUP_TABLE " << "default"<< endl;
        for (uint i = 0; i < cells.size() ; i++)
        {
            for ( uint j=0; j < cells.at(i).meshes.size() ; j++)
            {
                MeshOut << cells.at(i).meshes.at(j).concentrations2.at(2)<<" " <<cells.at(i).meshes.at(j).concentrations2.at(2) << " "<<cells.at(i).meshes.at(j).concentrations2.at(2) <<endl ;
            }
        }

        MeshOut << "SCALARS CLV3 " << "float"<< endl;
        MeshOut << "LOOKUP_TABLE " << "default"<< endl;
        for (uint i = 0; i < cells.size() ; i++)
        {
            for ( uint j=0; j < cells.at(i).meshes.size() ; j++)
            {
                MeshOut << cells.at(i).meshes.at(j).concentrations2.at(3)<<" " <<cells.at(i).meshes.at(j).concentrations2.at(3) << " "<<cells.at(i).meshes.at(j).concentrations2.at(3) <<endl ;
            }
        }
        
        MeshOut << "SCALARS ck " << "float"<< endl;
        MeshOut << "LOOKUP_TABLE " << "default"<< endl;
        for (uint i = 0; i < cells.size() ; i++)
        {
            for ( uint j=0; j < cells.at(i).meshes.size() ; j++)
            {
                MeshOut << cells.at(i).meshes.at(j).concentrations2.at(4)<<" " <<cells.at(i).meshes.at(j).concentrations2.at(4) << " "<<cells.at(i).meshes.at(j).concentrations2.at(4) <<endl ;
            }
        }

        MeshOut << "SCALARS ck_R " << "float"<< endl;
        MeshOut << "LOOKUP_TABLE " << "default"<< endl;
        for (uint i = 0; i < cells.size() ; i++)
        {
            for ( uint j=0; j < cells.at(i).meshes.size() ; j++)
            {
                MeshOut << cells.at(i).meshes.at(j).concentrations2.at(5)<<" " <<cells.at(i).meshes.at(j).concentrations2.at(5) << " "<<cells.at(i).meshes.at(j).concentrations2.at(5) <<endl ;
            }
        }
        

        MeshOut << "SCALARS CK " << "float"<< endl;
        MeshOut << "LOOKUP_TABLE " << "default"<< endl;
        for (uint i = 0; i < cells.size() ; i++)
        {
            for ( uint j=0; j < cells.at(i).meshes.size() ; j++)
            {
                MeshOut << cells.at(i).meshes.at(j).concentrations2.at(6)<<" " <<cells.at(i).meshes.at(j).concentrations2.at(6) << " "<<cells.at(i).meshes.at(j).concentrations2.at(6) <<endl ;
            }
        }
    }
    MeshOut.close();
}
 //---------------------------------------------------------------------------------------------
void SignalTissue::FullModel_Diffusion()
{
    for (unsigned int i=0; i<cells.size(); i++)
    {
        for (unsigned int j=0; j< cells.at(i).meshes.size(); j++)
        {
            cells.at(i).meshes.at(j).Flux.clear() ;
        }
    }
    
    for (unsigned int i =0 ; i<cells.size(); i++)
    {
        
        cells.at(i).FullModel_SelfDiffusion (cellType) ;
    }
    if (cellType == plant )
    {

        for (unsigned int i =0 ; i<cells.size(); i++)
        {
            for (unsigned int j =0; j < cells.at(i).meshes.size(); j++)
            {
                int cellID = cells.at(i).meshes.at(j).connection.first ;
                int meshID = cells.at(i).meshes.at(j).connection.second ;
                //cellID== -1 means the mesh is a boundary mesh (internal connections only )
                if (cellID == -1)
                {
                    continue ;
                }
                else if (i < cellID )
                {
                    double tmpChannel = cells.at(i).meshes.at(j).channel.at(2) ;
                    double length = cells.at(i).meshes.at(j).length.at(2) ;
                    vector<double> deltaU ;
                    deltaU.resize(cells.at(i).meshes.at(j).concentrations.size() ) ;
                   transform( cells.at(cellID).meshes.at(meshID ).concentrations.begin()+1, cells.at(cellID).meshes.at(meshID ).concentrations.begin()+5,
                              cells.at(i).meshes.at(j).concentrations.begin()+1, deltaU.begin()+1, linearConfig(1.0 ,-1.0) );
                    vector<double> flux ;
                    flux.clear() ;
                    flux.resize(deltaU.size() ) ;
                  transform(deltaU.begin()+1, deltaU.begin()+5, flux.begin()+1, productNum(tmpChannel/(length * length ) ) ) ;
                    transform(flux.begin()+1, flux.begin()+5, cells.at(i).meshes.at(j).diffusions.begin()+1, flux.begin()+1 , productVec()  ) ;
                    cells.at(i).meshes.at(j).Flux.push_back(flux) ;
                    transform(flux.begin()+1, flux.begin()+5, flux.begin()+1, productNum(-1.0) ) ;
                    cells.at(cellID).meshes.at(meshID).Flux.push_back( flux ) ;
        
                }
            }
        }
    }
    
    else if (cellType == wingDisc)
    {
        for (unsigned int i =0 ; i<cells.size(); i++)
        {
            for (unsigned int j =0; j < cells.at(i).meshes.size(); j++)
            {
                int cellID = cells.at(i).meshes.at(j).connection.first ;
                int meshID = cells.at(i).meshes.at(j).connection.second ;
                //cellID== -1 means the mesh is a boundary mesh (internal connections only )
                if (cellID == -1)
                {
                    continue ;
                }
                else if (i < cellID )
                {
                    double tmpChannel = cells.at(i).meshes.at(j).channel.at(2) ;
                    double length = cells.at(i).meshes.at(j).length.at(2) ;
                    vector<double> deltaU ;
                    deltaU.resize(cells.at(i).meshes.at(j).concentrations.size() ) ;
                    transform( cells.at(cellID).meshes.at(meshID ).concentrations.begin()+4, cells.at(cellID).meshes.at(meshID ).concentrations.begin()+5,
                              cells.at(i).meshes.at(j).concentrations.begin()+4, deltaU.begin()+4, linearConfig(1.0 ,-1.0) );
                    double dif = cells.at(i).meshes.at(j).diffusions.at(4) ;
                    vector<double> flux ;
                    flux.clear() ;
                    flux.resize(deltaU.size() ) ;
                    transform(deltaU.begin()+4, deltaU.begin()+5, flux.begin()+4, productNum(tmpChannel/(length * length ) ) ) ;
                    transform(flux.begin()+4, flux.begin()+5 , flux.begin()+4 , productNum( dif)  ) ;
                    cells.at(i).meshes.at(j).Flux.push_back(flux) ;
                    transform(flux.begin()+4, flux.begin()+5, flux.begin()+4, productNum(-1.0) ) ;
                    cells.at(cellID).meshes.at(meshID).Flux.push_back( flux ) ;
                    
                }
            }
        }
    }
    
}

//---------------------------------------------------------------------------------------------

void SignalTissue::FullModel_AllCellProductions()         //call once, initialize the constants
{
    for (unsigned int i=0; i<cells.size(); i++)
    {
        if (cellType == false)      //plant
        {
            if (cells.at(i).layer >= thres_layer && abs(cells.at(i).centroid.at(0) ) < thres_Production )
            {
                // CLV3 layer 3+
                cells.at(i).productionW = 1.0 ;
                cells.at(i).productionC = 1.0 ;
                // ck layer 4+
                cells.at(i).productionCk  = 1.0 ;
                cells.at(i).productionCkR  = 1.0 ;
                
            }
            else
            {
                cells.at(i).productionW = 0 ;
                cells.at(i).productionC = 0 ;
                cells.at(i).productionCk  = 0 ;
                cells.at(i).productionCkR  = 0 ;
            }
            cells.at(i).productionPMad  = 0 ;
        }
        else    //WingDisc
        {
            cells.at(i).productionW = 0 ;
            cells.at(i).productionC = 0 ;
            cells.at(i).productionCk  = 0.04 ;         //corresponds to Dpp
            cells.at(i).productionCkR  = 0.04 ;        //corresponds to Tkv
            cells.at(i).productionPMad  = 0.01 ;
        }
        cells.at(i).FullModel_ProductionCell() ;
    }
}
//---------------------------------------------------------------------------------------------
void SignalTissue::FullModelEulerMethod()
{
    double dt = cells.at(0).meshes.at(0).dt ;
    std:: stringstream stream ;
    stream <<fixed << setprecision(4) << dt  ;
    std::string dtString = stream.str();
    string txtFileName = "Histogram_"+ to_string(frameIndex) + "dt" + dtString + ".txt" ;
    ofstream histagram;
    histagram.open(txtFileName.c_str());
     
    bool state = false ;
    int l = 0 ;
    while (
           state==false
           && l<=300000
           )
    {
        double smallValue = 0.0005 ;
        FullModel_Diffusion() ;
    //    if (l%1000==0) cout<<l/1000<<endl ;
    //    if (l%100==0) ParaViewMesh(l/100) ;
        for (unsigned int i = 0; i < cells.size(); i++)
        {
            for (unsigned int j =0; j < cells.at(i).meshes.size(); j++)
            {
                cells.at(i).meshes.at(j).FullModel_Euler(cellType , TissueRadius,tissueCenter ) ;
            }
        }
        
        state = true ;
        for (unsigned int i = 0; i < cells.size(); i++)
        {
            for (unsigned int j =0; j < cells.at(i).meshes.size(); j++)
            {
                if (cellType == plant)
                {
                    if(
                       abs(cells.at(i).meshes.at(j).concentrations2.at(1) - cells.at(i).meshes.at(j).concentrations.at(1) )/(cells.at(i).meshes.at(j).concentrations.at(1) + smallValue) > smallValue * dt
                     ||  abs(cells.at(i).meshes.at(j).concentrations2.at(3) - cells.at(i).meshes.at(j).concentrations.at(3) )/(cells.at(i).meshes.at(j).concentrations.at(3) + smallValue) > smallValue * dt
                        || abs(cells.at(i).meshes.at(j).concentrations2.at(6) - cells.at(i).meshes.at(j).concentrations.at(6) )/(cells.at(i).meshes.at(j).concentrations.at(6) + smallValue) > smallValue * dt
                       )
                    {
                        state = false ;
                        break ;
                    }
                }
                else if (cellType== wingDisc)
                {
                    //condition is not finilized
                    if (abs(cells.at(i).meshes.at(j).concentrations2.at(4) - cells.at(i).meshes.at(j).concentrations.at(4) )/(cells.at(i).meshes.at(j).concentrations.at(4) + smallValue) > smallValue * dt   ||
                        abs(cells.at(i).meshes.at(j).concentrations2.at(7) - cells.at(i).meshes.at(j).concentrations.at(7) )/(cells.at(i).meshes.at(j).concentrations.at(7) + smallValue) > smallValue * dt
                        )
                    {
                        state = false ;
                        break ;
                    }
                }
                
            }
            if (state == false) break ;
        }
        for (unsigned int i = 0; i < cells.size(); i++)
        {
            for (unsigned int j =0; j < cells.at(i).meshes.size(); j++)
            {
                cells.at(i).meshes.at(j).UpdateU() ;
            }
        }
        if (l%10==0 && l<200)
        {
            Cal_AllCellConcentration() ;
            double overallC = accumulate(tissueLevelU.begin(), tissueLevelU.end(), 0.0) ;
            histagram<<l * dt<<'\t'<< overallC <<endl ;
        }
        
        if (l%100==0 && l>= 200 )
        {
            Cal_AllCellConcentration() ;
            double overallC = accumulate(tissueLevelU.begin(), tissueLevelU.end(), 0.0) ;
            histagram<<l * dt <<'\t'<< overallC <<endl ;
        }
        
        l++ ;
    }
    ParaViewMesh(frameIndex) ;
    cout<<"l is equal to "<<l << endl ;
    double value = 0 ;
    for (unsigned int j =0 ; j < cells.size();j++)
    {
        for (unsigned int i =0; i< cells.at(j).meshes.size(); i++)
        {
            value += cells.at(j).meshes.at(i).concentrations.at(4) ;
        }
    }
    cout<< "value is "<<value<<endl ;
  //  cout<< "number of steps needed is " << l <<endl ;
     
}
//---------------------------------------------------------------------------------------------

void SignalTissue::Cal_AllCellConcentration()
{
    tissueLevelConcentration.clear() ;
    tissueLevelU.clear() ;
    for (unsigned int i =0; i<cells.size(); i++)
    {
        cells.at(i).CellLevelConcentration(cellType) ;
       // cells.at(i).CellLevelConcentration2(cellType) ;
        if (equationsType== simpleODE)
        {
            tissueLevelU.push_back(cells.at(i).cellU ) ;
        }
        else
        {
            tissueLevelConcentration.push_back(cells.at(i).cellConcentration ) ;
            tissueLevelU.push_back(cells.at(i).cellConcentration.at(0) ) ;
        }
    }
    
}
//---------------------------------------------------------------------------------------------
void SignalTissue::Cal_ReturnSignal()
{
    if (equationsType == fullModel)
    {
        int sgnl ;
        if (cellType== wingDisc)
        {
            sgnl = 0 ;     //Dpp
        }
        else
        {
            sgnl = 3 ;     //CLV3
        }
        tissueLevelU.clear() ;
        for (unsigned int i=0 ; i < cells.size() ; i++)
        {
            tissueLevelU.push_back( tissueLevelConcentration.at(i).at(sgnl) ) ;
        }
    }
}
//---------------------------------------------------------------------------------------------
void SignalTissue:: Initialize_Concentrations(vector<vector<double> > oldConcentrations )
{
    for (unsigned int i=0; i<oldConcentrations.size() ; i++)
    {
        for (unsigned int j=0; j< oldConcentrations.at(i).size() ; j++)
        {
             if( isfinite(oldConcentrations.at(i).at(j) )== false  )
                {
                    return ;
                }
        }
    }
    
    // for now, concentrations have nothing to do with volume
    for (unsigned int i=0; i< cells.size(); i++)
    {
        int meshSize = static_cast<int>( cells.at(i).meshes.size() ) ;
        for (unsigned int j=0; j< cells.at(i).meshes.size(); j++)
        {
            if (cellType== plant)
            {
                transform(cells.at(i).meshes.at(j).concentrations.begin(), cells.at(i).meshes.at(j).concentrations.end()-1,
                          oldConcentrations.at(i).begin() , cells.at(i).meshes.at(j).concentrations.begin(), linearConfig(0.0,1.0) ) ;
                transform(cells.at(i).meshes.at(j).concentrations2.begin(), cells.at(i).meshes.at(j).concentrations2.end()-1,
                          oldConcentrations.at(i).begin() , cells.at(i).meshes.at(j).concentrations2.begin(), linearConfig(0.0,1.0) ) ;
            }
            // cellType == wingDisc
            else
            {
                transform(cells.at(i).meshes.at(j).concentrations.begin()+4, cells.at(i).meshes.at(j).concentrations.end() ,
                          oldConcentrations.at(i).begin() , cells.at(i).meshes.at(j).concentrations.begin()+4 , linearConfig(0.0,1.0) ) ;
                transform(cells.at(i).meshes.at(j).concentrations2.begin()+4, cells.at(i).meshes.at(j).concentrations2.end() ,
                          oldConcentrations.at(i).begin() , cells.at(i).meshes.at(j).concentrations2.begin()+4 , linearConfig(0.0,1.0) ) ;
            }
        }
    }
}

//---------------------------------------------------------------------------------------------

void SignalTissue::WriteConcentrations(string timer)
{
    if( writeVtk == false )
    {
        cout<< "WriteConcentratons is off" << endl ;
        return ;
    }
    string number = to_string(frameIndex) ;
    ofstream concentrationsData (timer + "Concentrations_"+ number +".txt") ;
    concentrationsData << cells.size()<<endl ;
    for (unsigned int i=0; i<tissueLevelConcentration.size(); i++)
    {
        concentrationsData<<i ;
        for (unsigned int j=0; j<tissueLevelConcentration.at(i).size(); j++)
        {
            concentrationsData<<'\t'<<tissueLevelConcentration.at(i).at(j) ;
        }
        concentrationsData<<endl ;
    }
    concentrationsData.close() ;
}
//---------------------------------------------------------------------------------------------
void SignalTissue::ReadConcentrations()
{
    
    if (cellType==wingDisc)
    {
        int a ;
        double b,c,d,e ;
        vector<double> data ;
        string number = to_string(frameIndex) ;
        ifstream concentrationsData ("newConcentrations_"+ number +".txt") ;
        if (concentrationsData.is_open())
        {
            concentrationsData >> a ;
            if (a==cells.size())
            {
                
                cout<<"concentrationsData is open"<<endl ;
                while (concentrationsData >> a >> b >> c >> d >> e)
                {
                    data = {b, c, d, e } ;
                    tissueLevelConcentration.push_back(data) ;
                }
                Initialize_Concentrations(tissueLevelConcentration) ;
            }
        }
            
    }
    else
    {
        int a ;
        double WmRNA,Wc,Wn,clv3,ck,ckr,CKR ;
        vector<double> data ;
        string number = to_string(frameIndex) ;
        ifstream concentrationsData ("newConcentrations_"+ number +".txt") ;
        if (concentrationsData.is_open())
        {
            concentrationsData >> a ;
            if (a==cells.size())
            {
                
                cout<<"concentrationsData is open"<<endl ;
                while (concentrationsData >> a >> WmRNA >> Wc >> Wn >> clv3 >> ck >> ckr >> CKR )
                {
                    data = { WmRNA , Wc, Wn, clv3, ck, ckr, CKR } ;
                    tissueLevelConcentration.push_back(data) ;
                }
                Initialize_Concentrations(tissueLevelConcentration) ;
            }
        }
    }
    
    
}
//---------------------------------------------------------------------------------------------

void SignalTissue::UpdateNanStatus()
{
    for (unsigned int i=0; i<tissueLevelConcentration.size() ; i++)
    {
        for (unsigned int j=0; j< tissueLevelConcentration.at(i).size() ; j++)
        {
            if( isfinite(tissueLevelConcentration.at(i).at(j) )== false  )
            {
                frameIsNan = true ;
                return ;
            }
        }
    }
}
//---------------------------------------------------------------------------------------------
void SignalTissue::WriteSignalingProfile()
{
    string number = to_string(frameIndex) ;
    std:: stringstream stream ;
    stream <<fixed << setprecision(1) << cells.at(0).meshes.at(0).diffusions.at(4)  ;
    std::string dif = stream.str();
    stream.str(string() ) ;
   // string dif = to_string(cells.at(0).meshes.at(0).diffusions.at(4) ) ;
    stream <<fixed << setprecision(3) << cells.at(0).meshes.at(0).degradations.at(4)  ;
    std::string deg = stream.str();
    stream.str(string() ) ;
  //  string deg = to_string(cells.at(0).meshes.at(0).degradations.at(4) ) ;
    stream <<fixed << setprecision(3) << cells.at(0).meshes.at(0).productions.at(4) ;
   // string pro = to_string(cells.at(0).meshes.at(0).productions.at(4) ) ;
    std::string pro = stream.str() ;
    stream.str(string() ) ;
    stream <<fixed << setprecision(5) << cells.at(0).meshes.at(0).dt ;
    std::string timeStep = stream.str() ;
    
    string cellTypeString = (cellType) ? "WingDisc" : "Plant" ;
    
    ofstream profile (cellTypeString + number + "D" + dif + "d" + deg + "p" + pro + "dt"+ timeStep + ".txt") ;
    
    for (int i=0 ; i< tissueLevelConcentration.size() ; i++)
    {
        profile << abs( cells.at(i).centroid.at(0)- tissueCenter.at(0) ) << '\t'
                << abs( cells.at(i).centroid.at(0)- tissueCenter.at(0) )/ TissueRadius ;
        for (unsigned int j = 0; j< tissueLevelConcentration.at(i).size() ; j++)
        {
            profile << '\t' << tissueLevelConcentration.at(i).at(j) ;
        }
        profile << endl ;
    }
    
    profile.close() ;
}

void SignalTissue::AddNoiseToChemical()
{

    for (unsigned int k=0 ;  k< cells.size(); k++)
    {
        double distYAbs = abs( cells.at(k).centroid.at(1) - tissueCenter.at(1) ) ;
        
        double dummy = (static_cast<double>(rand()) / RAND_MAX);
        double ranNum = NormalCDFInverse2(dummy);
        
        tissueLevelConcentration.at(k).at(0) *= 1.0 + (0.1*sin(0.2*3.141592*distYAbs)+0.12*ranNum) ;
    }

}
