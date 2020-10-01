#include "SignalCell.hpp"





//---------------------------------------------------------------------------------------------
SignalCell::SignalCell ()
{
    centroid.resize(2) ;
};
//---------------------------------------------------------------------------------------------

vector <vector<double> > SignalCell::Cal_NodeToNodeDist(vector <double> nodesXNeighbor, vector<double> nodesYNeighbor)
{
    vector<vector<double> > distNodeToNeighborNodes  ;
    for (unsigned int k=0; k < nodesX.size() ; k++)
    {
        vector <double> distNodeToNode ;
        distNodeToNode.clear() ;
        for (unsigned int l = 0; l < nodesXNeighbor.size(); l++)
        {
            distNodeToNode.push_back( Dist2D(nodesX.at(k),nodesY.at(k),nodesXNeighbor.at(l),nodesYNeighbor.at(l))) ;
        }
        
        distNodeToNeighborNodes.push_back(distNodeToNode) ;
    }
    return distNodeToNeighborNodes ;
}
//---------------------------------------------------------------------------------------------

void SignalCell::Find_NghbrCandidate ()
{
    double smallValue= 0.00001 ;
    for (unsigned int j = 0 ; j < cntrToCntr.size() ; j++)
    {
        if (cntrToCntr.at(j) < searchAreaForNghbr  && cntrToCntr.at(j) > smallValue)
        {
            nghbrCandidate.push_back(j) ;
        }
    }
}


//---------------------------------------------------------------------------------------------
void SignalCell::Cal_Centroid()
{
    
    centroid.at(0)= accumulate(nodesX.begin(), nodesX.end(), 0.0)/nodesX.size () ;
    centroid.at(1)= accumulate(nodesY.begin(), nodesY.end(), 0.0)/nodesY.size () ;
    
    
}

//---------------------------------------------------------------------------------------------
void SignalCell::Find_NghbrProperties ()
{
    for (unsigned int j = 0 ; j < nghbrCandidate.size(); j++)
    {
        vector<double> minDistnode ;
        vector<int> indexminDistnode ;
        double minDistCell ;
        //for all the nodes in the primary cell
        for (unsigned int k = 0; k < nodeDistoNghbrsCandidate.at(j).size(); k++)
        {
            minDistnode.push_back( *min_element(nodeDistoNghbrsCandidate.at(j).at(k).begin(),nodeDistoNghbrsCandidate.at(j).at(k).end() ) );
            
            indexminDistnode.push_back(static_cast<int>( min_element(nodeDistoNghbrsCandidate.at(j).at(k).begin(),nodeDistoNghbrsCandidate.at(j).at(k).end() )- nodeDistoNghbrsCandidate.at(j).at(k).begin() ) );
            
        }
        minDistCell = *min_element( minDistnode.begin(), minDistnode.end() ) ;
        if (minDistCell < thres_lateral)
        {
            Neighbor newNeighbor ;
            newNeighbor.CellID_Neighbor = nghbrCandidate.at(j) ;
            newNeighbor.minDistNode = minDistnode ;
            newNeighbor.nodeDistoNghbr = nodeDistoNghbrsCandidate.at(j) ;
            newNeighbor.indexMinDistNode = indexminDistnode  ;
            for (unsigned int l = 0; l < newNeighbor.minDistNode.size(); l++)
            {
                if (newNeighbor.minDistNode.at(l) < thres_lateral)
                {
                    newNeighbor.indexCellNode.push_back(l) ;
                    newNeighbor.indexNeighborNode.push_back(indexminDistnode.at(l)) ; // common nodes in the vector, can not trust on it
                }
            }
            neighbors.push_back(newNeighbor) ;
        }
    }
    
}
//---------------------------------------------------------------------------------------------
void SignalCell::NewEdge()
{
    vector<int> removeNode ;
    vector<double> addNodeX ;
    vector<double> addNodeY ;
    
    for (unsigned int j=0; j< neighbors.size(); j++)
    {
        removeNode.insert( removeNode.end(), neighbors.at(j).indexCellNode.begin() ,neighbors.at(j).indexCellNode.end() )  ;
        addNodeX.insert( addNodeX.end(), neighbors.at(j).intfX.begin(), neighbors.at(j).intfX.end()) ;
        addNodeY.insert( addNodeY.end(), neighbors.at(j).intfY.begin(), neighbors.at(j).intfY.end()) ;
    }
    sort(removeNode.begin(), removeNode.end() ) ;
    removeNode.erase( unique(removeNode.begin(),removeNode.end()) , removeNode.end() ) ;
    for (long k = removeNode.size()-1 ; k >= 0 ; k--)
    {
        nodesXNew.erase( nodesXNew.begin() + removeNode.at(k) ) ;
        nodesYNew.erase( nodesYNew.begin() + removeNode.at(k) ) ;
        noNeighboringNodesX.erase(noNeighboringNodesX.begin() + removeNode.at(k)) ;
        noNeighboringNodesY.erase(noNeighboringNodesY.begin() + removeNode.at(k)) ;
    }
    nodesXNew.insert( nodesXNew.end() , addNodeX.begin() , addNodeX.end() ) ;
    nodesYNew.insert( nodesYNew.end() , addNodeY.begin() , addNodeY.end() ) ;
}
//---------------------------------------------------------------------------------------------

void SignalCell::SortCCW () // vector<double> &verticesX , vector<double> &verticesY )
{
    vector<double> angle ;
    for (unsigned int j = 0; j < verticesX.size() ; j++)
    {
        angle.push_back( atan2(verticesY.at(j)-centroid.at(1), verticesX.at(j)-centroid.at(0)) ) ;
        angle.back() = fmod((angle.back() + 2* pi ) , (2*pi)) ;
    
    }
    vector<int > idx(angle.size());
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(),
         [&angle](int i1, int i2) {return angle[i1] < angle[i2];});
    vector<double> tmpX ;
    vector<double> tmpY ;
    for (unsigned int j = 0; j < verticesX.size(); j++)
    {
        tmpX.push_back(verticesX.at(idx.at(j) ) ) ;
        tmpY.push_back(verticesY.at(idx.at(j) ) ) ;
    }
    verticesX = tmpX ;
    verticesY = tmpY ;
}
//---------------------------------------------------------------------------------------------
void SignalCell::Cal_Vertices()
{
    for (unsigned int j =0; j < neighbors.size(); j++)
    {   int offset = 0 ;
        for (unsigned int k =0 ; k< neighbors.at(j).realIntersect.size(); k++)
        {
            if (neighbors.at(j).realIntersect.at(k) == true)
            {
                if (neighbors.at(j).CellID_Neighbor < neighbors.at(j).commonNeighbors.at(k).cellIDCommonNeighbor )
                {
                    verticesX.push_back( neighbors.at(j).intersectX.at(k - offset) ) ;
                    verticesY.push_back( neighbors.at(j).intersectY.at(k - offset) ) ;
                }
            }
            else
            {
                offset += 1 ;
            }
        }
    }
}
//---------------------------------------------------------------------------------------------
void SignalCell::Cal_Connections()
{
    for (unsigned int j=0 ; j< verticesX.size(); j++)
    {
        connections.push_back( (j+1) % verticesX.size() ) ;
    }
}

//---------------------------------------------------------------------------------------------
void SignalCell::Add_BoundaryVertice()
{
    if (boundary && noNeighboringNodesX.size() > 6)
    {
        double centerOfMassX = accumulate(noNeighboringNodesX.begin(),
                                          noNeighboringNodesX.end(), 0.0)/noNeighboringNodesX.size() ;
        double centerOfMassY =  accumulate(noNeighboringNodesY.begin(),
                                           noNeighboringNodesY.end(), 0.0)/noNeighboringNodesY.size() ;
        
        vector<double> tmpDist = Dist_PointToVec2D(centerOfMassX, centerOfMassY, noNeighboringNodesX, noNeighboringNodesY) ;
        int indexNode = min_element(tmpDist.begin(), tmpDist.end()) - tmpDist.begin() ;
        verticesX.push_back(noNeighboringNodesX.at(indexNode)) ;
        verticesY.push_back(noNeighboringNodesY.at(indexNode)) ;
    }
}
//---------------------------------------------------------------------------------------------
void SignalCell::Add_BoundaryVertice2()
{
    int thres = 15;
    if (noNeighboringNodesX.size() > thres)
    {
        vector<double> angle ;
        for (unsigned int j = 0; j < noNeighboringNodesX.size() ; j++)
        {
            angle.push_back( atan2(noNeighboringNodesY.at(j)-centroid.at(1), noNeighboringNodesX.at(j)-centroid.at(0)) ) ;
            angle.back() = fmod((angle.back() + 2*pi ), (2*pi)) ;
            
        }
        vector<int > idx(angle.size());
        iota(idx.begin(), idx.end(), 0);
        sort(idx.begin(), idx.end(),
             [&angle](int i1, int i2) {return angle[i1] < angle[i2];});
        vector<double> tmpX ;
        vector<double> tmpY ;
        for (unsigned int j = 0; j < noNeighboringNodesX.size(); j++)
        {
            tmpX.push_back(noNeighboringNodesX.at(idx.at(j) ) ) ;
            tmpY.push_back(noNeighboringNodesY.at(idx.at(j) ) ) ;
        }
        noNeighboringNodesX = tmpX ;
        noNeighboringNodesY = tmpY ;
        int number = ( noNeighboringNodesX.size()) / thres  ;
  //      cout<< noNeighboringNodesX.size() << '\t'<< thres << '\t' << number << endl ;
        int residue = (noNeighboringNodesX.size()) % thres  ;
        
        for (int j = 0; j < number ; j++)
        {
            verticesX.push_back(noNeighboringNodesX.at(thres*j + residue/2   )) ;
            verticesY.push_back(noNeighboringNodesY.at(thres*j + residue/2  )) ;
        }
    }
    
    
}
//---------------------------------------------------------------------------------------------
void SignalCell::Add_BoundaryVertice3 ()
{
    int thres = 5 ;
    vector<double> angle ;
    if (noNeighboringNodesX.size() > thres)
    {
        for (unsigned int j = 0; j < noNeighboringNodesX.size() ; j++)
        {
            angle.push_back( atan2(noNeighboringNodesY.at(j)-centroid.at(1), noNeighboringNodesX.at(j)-centroid.at(0)) ) ;
            angle.back() = fmod((angle.back()+ 2* pi ) , (2*pi)) ;
            
        }
        double preferedAngle = pi / 3 ;
        vector<vector<double> > PairwiseDist= Dist_VecToVec2D(noNeighboringNodesX, noNeighboringNodesY, noNeighboringNodesX, noNeighboringNodesY) ;
        vector<int > indices = Indices_MaxMatrix(PairwiseDist) ;
        pair<double , double > point ;  //(x,y)
        point.first = noNeighboringNodesX.at(indices.at(0)) ;
        point.second = noNeighboringNodesY.at(indices.at(0)) ;
        vector<double> tmpDist = Dist_PointToVec2D(point.first, point.second, verticesX, verticesY) ;
        int vertice1 = static_cast<int> (min_element(tmpDist.begin(), tmpDist.end()) - tmpDist.begin() ) ;
        point.first = noNeighboringNodesX.at(indices.at(1)) ;
        point.second = noNeighboringNodesY.at(indices.at(1)) ;
        tmpDist = Dist_PointToVec2D(point.first , point.second , verticesX, verticesY) ;
        int vertice2 = static_cast<int> (min_element(tmpDist.begin(), tmpDist.end()) - tmpDist.begin() ) ;
        double tetta1 = fmod(atan2(verticesY.at(vertice1) -centroid.at(1), verticesX.at(vertice1)-centroid.at(0)) + 2*pi ,2*pi) ;
        double tetta2 = fmod(atan2(verticesY.at(vertice2) -centroid.at(1), verticesX.at(vertice2)-centroid.at(0)) + 2*pi , 2*pi) ;
      //  double deltaTetta = fmod(tetta2 - tetta1 +2*pi, 2*pi) ;
        double deltaTetta = tetta2 - tetta1 ;

        // we assume the boundary starts from tetta1 and goes to tetta2 continuously
        if ((angle.at(0) < tetta1 && angle.at(0) > tetta2 ))
        {  //it is from tetta2 to tetta1
            double tmptetta = tetta1 ;
            tetta1 = tetta2 ;
            tetta2 = tmptetta ;
            deltaTetta = tetta2 - tetta1 ;
        }
        else if ( (angle.at(0) > tetta1 && angle.at(0) > tetta2) || (angle.at(0) < tetta1 && angle.at(0) < tetta2 ) )
        {   // periodic bounday condition
            if (deltaTetta > 0)
            {   //from tetta2 to tetta1
                double tmptetta = tetta1 ;
                tetta1 = tetta2 ;
                tetta2 = tmptetta ;
                deltaTetta = 2*pi + tetta2 - tetta1 ;
            }
            else
            {
                deltaTetta += 2*pi ;
            }
        }
        if (deltaTetta < 0.0)
        {
            cout<< "OMG!!! deltaTetta is negative :"<<deltaTetta<<endl ;
        }
        int number = round(deltaTetta / preferedAngle);
        double step = deltaTetta / number ;
        for (int i =1; i< number; i++)
        {
            double tmpDegree =fmod(tetta1 +i* step + 2*pi,2*pi ) ;
            vector<double> dtetta = Dist_pointToVec1D(tmpDegree, angle) ;
            int id = static_cast<int>( min_element(dtetta.begin(), dtetta.end()) - dtetta.begin() ) ;
            verticesX.push_back(noNeighboringNodesX.at(id)) ;
            verticesY.push_back(noNeighboringNodesY.at(id)) ;
        }
    }
    return ;
}

//---------------------------------------------------------------------------------------------
void SignalCell::Refine_NoBoundary ()
{
    
    vector<double> tmp ;
    vector<int> removeNoBoundaryNode ;
    int size = static_cast<int>(noNeighboringNodesX.size()) ;
    if (verticesX.size() > 0)
    {
        for (int i = 0 ;i < size ; i++)
        {
            tmp = Dist_PointToVec2D(noNeighboringNodesX.at(i), noNeighboringNodesY.at(i), verticesX, verticesY) ;
            double tmpMin = *min_element(tmp.begin(), tmp.end()) ;
            if (tmpMin < thres_refineNodes )
            {
                removeNoBoundaryNode.push_back(i) ;
            }
            
        }
        sort(removeNoBoundaryNode.begin(), removeNoBoundaryNode.end()) ;
        removeNoBoundaryNode.erase(unique(removeNoBoundaryNode.begin(), removeNoBoundaryNode.end()), removeNoBoundaryNode.end()) ;
        for (long j= removeNoBoundaryNode.size()-1 ; j>= 0 ; j--)
        {
            noNeighboringNodesX.erase(noNeighboringNodesX.begin() + removeNoBoundaryNode.at(j) ) ;
            noNeighboringNodesY.erase(noNeighboringNodesY.begin() + removeNoBoundaryNode.at(j) ) ;
            
        }
    }
}
//---------------------------------------------------------------------------------------------
void SignalCell::Refine_NodeXNew ()
{
    vector<double> tmp ;
    vector<int> removeNodeXNew ;
    int size = static_cast<int>(nodesXNew.size()) ;
    if (verticesX.size() > 0)
    {
        for (unsigned int i = 0 ;i < size ; i++)
        {
            tmp = Dist_PointToVec2D(nodesXNew.at(i), nodesYNew.at(i), verticesX, verticesY) ;
            double tmpMin = *min_element(tmp.begin(), tmp.end()) ;
            if (tmpMin < thres_refineNodes )
            {
                removeNodeXNew.push_back(i) ;
            }
            
        }
        sort(removeNodeXNew.begin(), removeNodeXNew.end()) ;
        removeNodeXNew.erase(unique(removeNodeXNew.begin(), removeNodeXNew.end()), removeNodeXNew.end()) ;
        for (long j=removeNodeXNew.size()-1 ; j>= 0 ; j--)
        {
            nodesXNew.erase(nodesXNew.begin() + removeNodeXNew.at(j) ) ;
            nodesYNew.erase(nodesYNew.begin() + removeNodeXNew.at(j) ) ;
            
        }
    }
}

//---------------------------------------Meshes------------------------------------------------
//---------------------------------------------------------------------------------------------
void SignalCell::Find_Mesh()
{
    int size = static_cast<int>(verticesX.size()) ;
    for (unsigned int i = 0; i < verticesX.size(); i++)
    {
        SignalMesh tmpMesh ;
        tmpMesh.triangleX.push_back( centroid.at(0) ) ;
        tmpMesh.triangleY.push_back( centroid.at(1) ) ;
        
        tmpMesh.triangleX.push_back(verticesX.at(i) ) ;
        tmpMesh.triangleY.push_back(verticesY.at(i) ) ;
        
        tmpMesh.triangleX.push_back(verticesX.at( (i+1) % size )) ;
        tmpMesh.triangleY.push_back(verticesY.at( (i+1) % size )) ;
        
        tmpMesh.Cal_MeshCenter() ;
        meshes.push_back(tmpMesh) ;
        
    }
    for (unsigned int i =0; i< meshes.size(); i++)
    {
        int size = static_cast<int>(meshes.size()) ;
        double tmpChannel = Dist2D(meshes.at(i).triangleX.at(0), meshes.at(i).triangleY.at(0),
                             meshes.at(i).triangleX.at(2), meshes.at(i).triangleY.at(2)) ;
        double length = Dist2D(meshes.at(i).center.at(0), meshes.at(i).center.at(1),
                               meshes.at((i+1) % size ).center.at(0), meshes.at((i+1) % size).center.at(1)) ;
        meshes.at(i).channel.at(0) = tmpChannel ;
        meshes.at((i+1) % size).channel.at(1) = tmpChannel ;
        meshes.at(i).length.at(0)= length ;
        meshes.at((i+1) % size).length.at(1) = length ;
    }
    for (unsigned int i = 0; i < meshes.size(); i++)
    {
        double tmpChannel = Dist2D(meshes.at(i).triangleX.at(1), meshes.at(i).triangleY.at(1),
                             meshes.at(i).triangleX.at(2), meshes.at(i).triangleY.at(2)) ;
        meshes.at(i).channel.at(2) = tmpChannel ;
    }
}

//---------------------------------------------------------------------------------------------
void SignalCell::Self_Diffusion()
{
    double dSelf = 4.0 ;
    for (unsigned int i =0; i< meshes.size(); i++)
    {
        int size = static_cast<int>(meshes.size()) ;
        double tmpChannel = meshes.at(i).channel.at(0) ;
        double length = meshes.at(i).length.at(0) ;
        double deltaU = meshes.at((i+1) % size).u1 - meshes.at(i).u1 ;
        double flux = dSelf * tmpChannel * deltaU / (length * length)  ;
        meshes.at(i).flux.push_back( flux ) ;
        meshes.at((i+1) % size).flux.push_back(- flux ) ;
    }
}
//---------------------------------------------------------------------------------------------
void SignalCell::FullModel_SelfDiffusion (bool type)
{
    if (type == false) //plant
    {
        for (unsigned int i =0; i< meshes.size(); i++)
        {
            int size = static_cast<int>(meshes.size()) ;
            double tmpChannel = meshes.at(i).channel.at(0) ;
            double length = meshes.at(i).length.at(0) ;
            vector<double> deltaU ;
            deltaU.resize(meshes.at(i).concentrations.size() ) ;
            transform( meshes.at((i+1)% size ).concentrations.begin()+1, meshes.at((i+1)% size ).concentrations.begin()+5,
                      meshes.at(i).concentrations.begin()+1, deltaU.begin()+1, linearConfig(1.0 ,-1.0) );
            vector<double> flux ;
            flux.resize(deltaU.size() ) ;
            transform(deltaU.begin()+1, deltaU.begin()+5, flux.begin()+1, productNum(tmpChannel/ (length * length ) ) ) ;
            transform(flux.begin()+1, flux.begin()+5, meshes.at(i).selfDiffusions.begin()+1, flux.begin()+1, productVec() ) ;
            meshes.at(i).Flux.push_back(flux) ;
            transform(flux.begin()+1, flux.begin()+5, flux.begin()+1, productNum(-1.0) ) ;
            meshes.at((i+1) % size).Flux.push_back( flux ) ;
            
        }
    }
    else
    {
        for (unsigned int i =0; i< meshes.size(); i++)
        {
            int size = static_cast<int>(meshes.size()) ;
            double tmpChannel = meshes.at(i).channel.at(0) ;
            double length = meshes.at(i).length.at(0) ;
            vector<double> deltaU ;
            deltaU.resize(meshes.at(i).concentrations.size() ) ;
            transform( meshes.at((i+1)% size ).concentrations.begin()+4, meshes.at((i+1)% size ).concentrations.begin()+5,
                      meshes.at(i).concentrations.begin()+4, deltaU.begin()+4, linearConfig(1.0 ,-1.0) );
            vector<double> flux ;
            flux.resize(deltaU.size() ) ;
            transform(deltaU.begin()+4, deltaU.begin()+5, flux.begin()+4, productNum(tmpChannel/( length * length ) ) ) ;
            transform(flux.begin()+4, flux.begin()+5, meshes.at(i).selfDiffusions.begin()+4, flux.begin()+4, productVec() ) ;
            meshes.at(i).Flux.push_back(flux) ;
            transform(flux.begin()+4, flux.begin()+5, flux.begin()+4, productNum(-1.0) ) ;
            meshes.at((i+1) % size).Flux.push_back( flux ) ;
        }
    }
}
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void SignalCell:: FullModel_ProductionCell()
{
    for (unsigned int j=0; j< meshes.size(); j++)
    {
        meshes.at(j).productions = { productionW , 0.0, 0.0, productionC, productionCk, productionCkR, 0.0 , productionPMad } ;
    }
}
//---------------------------------------------------------------------------------------------
void SignalCell:: CellLevelConcentration(bool type)
{
    //Calculate average concentration based on area of each mesh
    cellConcentration.clear() ;
    cellU = 0.0 ;
    if (type == false)  //plant
    {
        cellConcentration.resize(meshes.at(0).concentrations.size() -1 ) ;
    }
    else
    {
        cellConcentration.resize(meshes.at(0).concentrations.size() -4  ) ;
    }
    
    for (unsigned int j=0 ; j< meshes.size() ; j++ )
    {
        cellU += meshes.at(j).u1 ;
        if (type == false)  //plant
        {
        transform(meshes.at(j).concentrations.begin(), meshes.at(j).concentrations.end()-1 , cellConcentration.begin(), cellConcentration.begin() , linearConfig (meshes.at(j).area,1.0) ) ;
        }
        else    //wingDisc
        {
            transform(meshes.at(j).concentrations.begin()+4, meshes.at(j).concentrations.end() , cellConcentration.begin(), cellConcentration.begin() , linearConfig (meshes.at(j).area, 1.0) ) ;
        }
    }
    transform(cellConcentration.begin(), cellConcentration.end(), cellConcentration.begin(), productNum(1.0/areaCell ) ) ;
}

//---------------------------------------------------------------------------------------------
void SignalCell:: CellLevelConcentration2(bool type)
{
    //Calculate average concentration based on number of the meshes on each cell
    cellConcentration.clear() ;
    cellU = 0.0 ;
    //  cellConcentration.resize(meshes.at(0).concentrations.size() ) ;
    if (type == false)  //plant
    {
        cellConcentration.resize(meshes.at(0).concentrations.size() -1 ) ;
    }
    else
    {
        cellConcentration.resize(meshes.at(0).concentrations.size() -4  ) ;
    }
    
    for (unsigned int j=0 ; j< meshes.size() ; j++ )
    {
        cellU += meshes.at(j).u1 ;
        if (type == false)  //plant
        {
            transform(meshes.at(j).concentrations.begin(), meshes.at(j).concentrations.end()-1 , cellConcentration.begin(), cellConcentration.begin() , linearConfig (1.0,1.0) ) ;
        }
        else    //wingDisc
        {
            transform(meshes.at(j).concentrations.begin()+4, meshes.at(j).concentrations.end() , cellConcentration.begin(), cellConcentration.begin() , linearConfig (1.0, 1.0) ) ;
           
        }
    }
    transform(cellConcentration.begin(), cellConcentration.end(), cellConcentration.begin(), productNum(1.0/meshes.size() ) ) ;
}
