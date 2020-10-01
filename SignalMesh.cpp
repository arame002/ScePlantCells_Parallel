
#include "SignalMesh.hpp"

SignalMesh::SignalMesh()
{
    connection.first = -1 ;
    connection.second = -1 ;
    channel.resize(3) ;
    length.resize(3) ;
    concentrations.resize(8) ;
    concentrations2.resize(8) ;
    diffusions.resize(8) ;
    degradations.resize(8) ;
    productions.resize(8) ;
    rates.resize(8) ;
    D = 20.0 ;
    d = 1.8 ;
    diffusions = {0,D,0,D,D,0,0,0} ;
    selfDiffusions = {0,D,0,D,D,0,0,0} ;
    degradations = {d, d, d, d, d, d, 2.0*d , d } ;
    concentrations = {0,0,0,0,0,0,0,0} ;
    concentrations2 = {0,0,0,0,0,0,0,0} ;
  //  productions = {1,0,0,1,1,1,0,1} ;
    rates = {1,1,1,0,0.003,0.003,1,0} ; // <r_c, r_im, r_ex, 0 , k_on , k_on, k_off>
    powers.resize(4) ;
    powers = {-2,1,2,1} ;
    dt = 0.01/ D ;
    
    
   // dpp.resize(4)  ;
   // dpp2.resize(4) ;
   // diff = 1.0 ;
   // degreDpp.resize(4)  ;
   // producDpp.resize(4) ;
}

//--------------------------------------------------------------------------------------------
void SignalMesh::Cal_MeshCenter()
{
    
    center.push_back(accumulate(triangleX.begin(), triangleX.end(), 0.0)/triangleX.size ()) ;
    center.push_back(accumulate(triangleY.begin(), triangleY.end(), 0.0)/triangleY.size ()) ;
    
}
//---------------------------------------------------------------------------------------------
double SignalMesh::Cal_MeshArea()
{
    pair<double, double> vec1 ;
    vec1.first = triangleX.at(1) - triangleX.at(0) ;
    vec1.second = triangleY.at(1) - triangleY.at(0) ;
    
    pair<double, double> vec2 ;
    vec2.first = triangleX.at(2) - triangleX.at(0) ;
    vec2.second = triangleY.at(2) - triangleY.at(0) ;
    area = TriangleArea(vec1.first, vec1.second, vec2.first, vec2.second) ;
    return area ;
    
}
//---------------------------------------------------------------------------------------------
void SignalMesh::Cal_UStar ()
{
    double totalFlux = accumulate(flux.begin(), flux.end(), 0.0) ;
    uStar = u1 + dt * (totalFlux + c - d*u1) ;
}
//---------------------------------------------------------------------------------------------
void SignalMesh::Euler()
{
    double totalFlux = accumulate(flux.begin(), flux.end(), 0.0) ;
    u2 = u1 + dt * (totalFlux + c - d * u1 ) ;
    flux.clear() ;
}
//---------------------------------------------------------------------------------------------
void SignalMesh::UpdateU()
{
    u1 = u2   ;
    concentrations = concentrations2 ;
}
//---------------------------------------------------------------------------------------------
void SignalMesh::FullModel_Euler(bool type, double radius, vector<double> tissueCntr)
{
    
    vector<double> degradationChanges = UpdateDegradation(type, radius,tissueCntr ) ;
    vector<double> rateChanges = UpdateRates(type) ;
    
    
    vector<double> totalFlux ;
    totalFlux.clear() ;
    
    for (unsigned int j= 0 ; j< diffusions.size() ; j++ )
    {
       totalFlux.push_back( sum_over_vec(Flux, j) ) ;
        
    }
    if (type == false)  //plant
    {
    
        // Diffusion changes from 1st to 4th
        transform(concentrations.begin()+1, concentrations.begin()+5, totalFlux.begin()+1, concentrations2.begin()+1, linearConfig (1.0, dt ) ) ;
        // degradation changes
        transform(concentrations2.begin(), concentrations2.end()-1, degradationChanges.begin() ,concentrations2.begin(), linearConfig( 1.0,-dt )  ) ;
        //rate changes
        transform(concentrations2.begin()+1 , concentrations2.end()-1, rateChanges.begin()+1 , concentrations2.begin()+1 , linearConfig(1.0 , -dt ) ) ;
        concentrations2.at(1) += (rateChanges.at(0) + rateChanges.at(2) ) * dt ;
        concentrations2.at(2) += rateChanges.at(1) * dt ;
        concentrations2.at(4) += rateChanges.at(6) * dt ;
        concentrations2.at(5) += rateChanges.at(6) * dt ;
        concentrations2.at(6) += rateChanges.at(4) * dt ;
        
        //production changes
        vector<double> tmpProduction = productions ;
        tmpProduction.at(0) *= 1/(1+ pow(concentrations.at(3)/kcw,powers.at(0) ) ) ;
        tmpProduction.at(3) *= UpdateCp() ;
        transform(concentrations2.begin(), concentrations2.end()-1, tmpProduction.begin(), concentrations2.begin(),linearConfig(1.0,dt) ) ;
    }
    else    // WingDisc
    {
        // Diffusion changes only 4th
        transform(concentrations.begin()+4, concentrations.begin()+5, totalFlux.begin()+4, concentrations2.begin()+4, linearConfig (1.0, dt ) ) ;
        // degradation changes (from 4th to the end(7th)
        transform(concentrations2.begin()+4, concentrations2.end(),degradationChanges.begin()+4 ,concentrations2.begin()+4, linearConfig( 1.0,-dt )  ) ;
        //rate changes
        /*  test
        transform(concentrations2.begin()+4 , concentrations2.end()-1, rateChanges.begin()+4 , concentrations2.begin()+4 , linearConfig(1.0 , - dt ) ) ;
        concentrations2.at(4) += rateChanges.at(6) * dt ;
        concentrations2.at(5) += rateChanges.at(6) * dt ;
        concentrations2.at(6) += rateChanges.at(4) * dt ;
        */
        //production changes
        vector<double> tmpProduction = productions ;
        double rs = 0.12 ;
        int ns = 40 ;
        double cntX = tissueCntr.at(0) ;
        double cntY = tissueCntr.at(1) ;
        tmpProduction.at(4) *= 1/(1+ pow(abs(triangleX.at(0)-cntX )/(rs * radius), ns ) ) ;        //check power
        tmpProduction.at(5) *= (cMax-cMin)/(1+ pow(concentrations.at(7)/kP,powers.at(2) ) ) ;        //check power
        tmpProduction.at(5) += cMin * productions.at(5) ;
        tmpProduction.at(5) *= 1/(1+ pow(abs(triangleX.at(0)-cntX )/(rs * radius) , -ns ) ) ;        //check power
        tmpProduction.at(5) += (3.0/6.0) * productions.at(5) ;
        tmpProduction.at(7) *= 1/(1+ pow(concentrations.at(6)/kLR , powers.at(0) ) ) ;        //check power
        transform(concentrations2.begin()+4 , concentrations2.end(), tmpProduction.begin()+4, concentrations2.begin()+4, linearConfig(1.0 , dt ) ) ;
        
    }
    
    
}
//---------------------------------------------------------------------------------------------
vector<double> SignalMesh::UpdateDegradation(bool type , double radius, vector<double> tissueCntr)
{
    int ns=40 ;
    double rs = 0.12 ;
    double cntX = tissueCntr.at(0) ;
    vector<double> tmp ;
    tmp.resize(degradations.size() ) ;
    if (type == false) //plant
    {
        
        degradations.at(1) = degradations.at(0) * (cMin + (cMax - cMin)/(1 + pow(concentrations.at(6)/ kCK,powers.at(1)) ) )
                                                * (eMin + (eMAx- eMin)/(1+pow(concentrations.at(3)/kcw,1)  ) ) ;
     //  degradations.at(2) = ( degradations.at(0)/(1 + pow(concentrations.at(6)/ kCK,powers.at(1) ) ) ) *
     //                                                   (1 / (1+pow(concentrations.at(3)/kcw, powers.at(3) ) ) ) ;
        degradations.at(2) = degradations.at(0) / (1 + pow(concentrations.at(2)/kcw,powers.at(3) ) ) ;
        transform(concentrations.begin(), concentrations.end()-1, degradations.begin(), tmp.begin(), productVec() ) ;
    }
    else    // Wing_Disc
    {
        tmp = degradations ;
        transform(concentrations.begin()+4, concentrations.end() , degradations.begin()+4, tmp.begin()+4, productVec() ) ;
        //test
     //   tmp.at(4) *= 1 / (1+pow(abs(triangleX.at(0)- cntX )/(rs * radius) , ns) ) ;
    }
    return tmp ;

}
//---------------------------------------------------------------------------------------------

vector<double> SignalMesh:: UpdateRates(bool type )
{
    vector<double> tmp ;
    tmp.resize(concentrations.size()) ;
    
    if (type == false) //plant
    {
        // <r_c, r_im, r_ex, 0 , k_on , k_on, k_off>
        transform(concentrations.begin(), concentrations.end()-1,rates.begin() , tmp.begin(), productVec() ) ;
        tmp.at(2) *= 1/(1 + pow(concentrations.at(3)/kcw, powers.at(2) ) ) ;
        tmp.at(2) *= 1/(1 + pow(concentrations.at(6)/kCK, powers.at(1) ) ) ;
        tmp.at(4) *= concentrations.at(5) ;
        tmp.at(5) = tmp.at(4) ;

    }
    else
    {
        transform(concentrations.begin()+4, concentrations.end()-1,rates.begin()+4 , tmp.begin()+4, productVec() ) ;
        tmp.at(4) *= concentrations.at(5) ;
        tmp.at(5) = tmp.at(4) ;
    }
    return tmp ;
}
//---------------------------------------------------------------------------------------------
double SignalMesh::UpdateCp ()
{
    double w1 = 10 ;
    double w2 = 100 ;
    double kLow = 1 ;
    double wn = concentrations.at(2) ;
    double cp = 1 ;
    if (wn > w1)
    {
        cp *= 2/(1+(wn-w1)/(w2-w1) ) ;
    }
    else
    {
        cp *= 2/(1+(w1-wn)/kLow) ;
    }
    return cp ;
}

//---------------------------------------------------------------------------------------------

