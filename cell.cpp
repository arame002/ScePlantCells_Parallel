//cell.cpp:
//===================
// Forward Declarations
// Include Dependencies
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <random>
#include <memory>
#include "phys.h"
#include "coord.h"
#include "node.h"
#include "rand.h"
#include "cell.h"
#include "tissue.h"
//===================

// Cell Class Member functions

// Constructors
// this constructor is used in
// the divsion function
Cell::Cell(Tissue* tissue, int parent) {
	my_tissue = tissue;
	//rank assigned in division function
	//layer inherited from parent	
	//boundary cells don't divide
	boundary = 0;
	ancestry = parent;
	//stem cells don't divide 
	stem = 0;
	//damping assigned in div function
	//just divided so reset life length
	life_length = 0;
	//cyt nodes reassigned in division function
	//start at zero
	num_cyt_nodes = 0;
	//wall nodes renumbered in division function
	//start at zero
	num_wall_nodes = 0;
	Cell_Progress = 0;
	//center calculate in division function
	//will calculate signals in div function
	wuschel = 0;
	cytokinin = 0;
	//growth_rate assigned in div function
	//growth direction assigned in division
	//neighbors assigned in div function
	//left corner assigned in division
}
//this constructor is used to initialize first set of cells
//calls set_growth_rate which detemrines growth rate based on WUS CONC
Cell::Cell(int rank, Coord center, double radius, Tissue* tiss, int layer, int boundary, int stem)    {
	this->my_tissue = tiss;
	this->rank = rank;
	this->ancestry = rank;
	this->layer = layer;
	//if boundary is equal to one 
	//then the cell will have higher damping
	//which is assigned below
	//boundary conditions are read in from initial text file
	this->boundary = boundary;
	this-> stem = stem;
	//set damping for cells that act as anchor points
	if(stem == 1) {
		this->damping = BOTTOM_DAMP;
	} else if((this->boundary == 1)) {
		this->damping =  BOUNDARY_DAMP;
	} else {
		this->damping = REG_DAMP;
	}
	life_length = 0;
	//cyt nodes initialized in tissue constructor which
	//calls the makes nodes function on each new cell
	num_cyt_nodes = 0;
	//wall nodes initialized in tissue constructor which 
	//calls the make nodes function on each new cell
	num_wall_nodes = 0;
	// REAL 
	//Cell_Progress = unifRandInt(0,10);
	// TEST - QUICK DIVISION
	Cell_Progress = unifRandInt(7,10);
	this->cell_center = center;
	this->calc_WUS();
	this->calc_CK();
	this->set_growth_rate();
	if(this->boundary == 1) {
		this->growth_direction = Coord(0,0);
	} else if(this->stem == 1) {
		this->growth_direction = Coord(0,1);
	} else if((this->layer == 1)||(this->layer == 2)) {
		//Growth Direction:  Apical growth in layer 1 is taken to be Isotropic
		//this->growth_direction = Coord(0,0);
		//Set to be biased for the algorithm.
		//Vertical
		this->growth_direction = Coord(0,1);
		//Horizontal
		//this->growth_direction = Coord(1,0);
	} else {
		this->growth_direction = Coord(0,0);
	}
	//cout << "layer" << this->layer << endl;
	//cout << "stem" << this->stem << endl;
	//cout << "boundary" << this-> boundary << endl;
	//cout << "gd" << this->growth_direction << endl;
	//cout << "damping" << this->damping << endl;
	//neighbors update function called after initialization
	//left corner set in make nodes function called by tissue constuctor
}
//calls compute membrane equi length for each node
//calls compute linear spring constant for each node
//calls compute bending spring constant for each node
//calls update wall equi angles for each node
//calls update wall angles to get initial angle of each node
void Cell::make_nodes(double radius) {

	//assemble the membrane
	int num_Init_Wall_Nodes = Init_Wall_Nodes;
	double angle_increment = (2*pi)/num_Init_Wall_Nodes;

	//make all wall nodes
	double curr_X;
	double curr_Y;
	Coord location = this->cell_center;;
	double curr_theta = 0;
	curr_X = cell_center.get_X() + radius*cos(curr_theta);
	curr_Y = cell_center.get_Y() + radius*sin(curr_theta);
	location = Coord(curr_X,curr_Y);
	//make the first node
	shared_ptr<Cell> this_cell = shared_from_this();
	shared_ptr<Wall_Node> prevW = make_shared<Wall_Node>(location,this_cell);
	shared_ptr<Wall_Node> currW = prevW;
	wall_nodes.push_back(prevW);
	num_wall_nodes++;
	shared_ptr<Wall_Node> orig(prevW);
	//this will be the "starter" node
	this->left_Corner = orig;

	//make successive nodes
	for(int i = 0; i<num_Init_Wall_Nodes-1; i++) {
		curr_theta = curr_theta + angle_increment;
		curr_X = cell_center.get_X() + radius*cos(curr_theta);
		curr_Y = cell_center.get_Y() + radius*sin(curr_theta);
		location = Coord(curr_X,curr_Y);
		shared_ptr<Wall_Node> new_node =make_shared<Wall_Node>(location,this_cell);
		currW = new_node;
		wall_nodes.push_back(currW);
		num_wall_nodes++;
		prevW->set_Left_Neighbor(currW);
		currW->set_Right_Neighbor(prevW);
		prevW = currW;
	}

	//connect last node to starter node
	currW->set_Left_Neighbor(orig);
	orig->set_Right_Neighbor(currW);

	//reindex wall node vector for vertical growing cells
	//this is old and  not necessary
	//if(this->growth_direction == Coord(0,1)){
	//	for(int i = 0; i<35; i++){
	//		this->left_Corner = left_Corner->get_Left_Neighbor();
	//	}
	//shared_ptr<Wall_Node> curr= left_Corner;
	//wall_nodes.clear();
	//do{
	//	wall_nodes.push_back(curr);
	//	curr = curr->get_Left_Neighbor();
	//} while(curr != left_Corner);
	//}

	//here is where most private member variables are set
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	//damping inherited from cell
	double new_damping = this->get_Damping();
	//variable to hold membr_equ_len
	double l_thresh = 0;
	double k_lin = 0;
	double k_bend = 0;
	for(unsigned int i = 0; i < walls.size();i++) {	
		walls.at(i)->set_Damping(new_damping);
		l_thresh = compute_membr_thresh(walls.at(i));
		k_lin = compute_k_lin(walls.at(i));
		k_bend = compute_k_bend(walls.at(i));
		walls.at(i)->set_membr_len(l_thresh);
		walls.at(i)->set_K_LINEAR(k_lin);
		walls.at(i)->set_K_BEND(k_bend);
	}

	//insert cytoplasm nodes
	int num_init_cyt_nodes = Init_Num_Cyt_Nodes + Cell_Progress;
	this->Cell_Progress = num_init_cyt_nodes;
	double scal_x_offset = 0.8;
	//Coord location;
	double x;
	double y;
	for (int i = 0; i < num_init_cyt_nodes; i++) {
		// USING POSITIONS OF CELL CENTER FOR CYT NODE ALLOCATION
		// ---distributes more evenly throughout start cell
		double rand_radius = (static_cast<double>(rand()) / RAND_MAX)*scal_x_offset*radius;
		double rand_angle = (static_cast<double>(rand()) / RAND_MAX)*2*pi;
		x = cell_center.get_X()+ rand_radius*cos(rand_angle);
		y = cell_center.get_Y()+ rand_radius*sin(rand_angle);
		location = Coord(x,y);

		shared_ptr<Cyt_Node> cyt = make_shared<Cyt_Node>(location,this_cell);
		cyt->set_Damping(new_damping);
		cyt_nodes.push_back(cyt);
		num_cyt_nodes++;
	}

	//update equilibrium angle
	update_Wall_Equi_Angles();
	//update wall angles
	update_Wall_Angles();	
	//is_divided = false;
	return;
}
// Destructor
Cell::~Cell() {
	//not needed using smartpointers
}
//=============================================================
//========================================
// Getters and Setters
//========================================
//=============================================================
void Cell::set_Rank(const int id) {
	this->rank = id;
	return;
}
void Cell::set_Layer(int layer) {
	this->layer = layer;
	return;
}
void Cell::set_Damping(double new_damping) {
	this->damping = new_damping;
	return;
}
void Cell::update_Life_Length() {
	life_length++;
	return;
}
void Cell::reset_Life_Length(){
	life_length = 0;
	return;
}
int Cell::get_Node_Count() {
	return num_wall_nodes + num_cyt_nodes;
}
void Cell::get_Wall_Nodes_Vec(vector<shared_ptr<Wall_Node>>& walls) {
	walls = this->wall_nodes;
	return;
}
void Cell::add_wall_node_vec(shared_ptr<Wall_Node> curr) {
	this->wall_nodes.push_back(curr);
	this->num_wall_nodes++;
	return;
}
void Cell::get_Cyt_Nodes_Vec(vector<shared_ptr<Cyt_Node>>& cyts) {
	cyts = cyt_nodes;
	return;
}
void Cell::update_cyt_node_vec(shared_ptr<Cyt_Node> new_node){
	this->cyt_nodes.push_back(new_node);
	this->num_cyt_nodes++;
	return;
}
void Cell::reset_Cell_Progress(){
	this->Cell_Progress = 0;
	return;
}
void Cell::update_Cell_Progress() {
	this->Cell_Progress++;
	return;
}
void Cell::calc_WUS() {
	this->wuschel = 109.6*exp(-0.02928*(cell_center - Coord(0,0)).length()) + 27.69*exp(-0.0008808*(cell_center - Coord(0,0)).length());
	//TEST - Reduce WUS since this will give painfully slow growth for one_cell. Divided by 4 to get center WUS around 32.
	this->wuschel = this->wuschel / static_cast<double>(4);
	return;
}
void Cell::calc_CK() {
	this->cytokinin = 132.9*exp(-0.01637*(cell_center-Coord(0,-40)).length());
	return;
}
void Cell::set_growth_rate() {
	// ULTRAFAST - TEST (Potentially unstable)
	//this->growth_rate = unifRandInt(2500,3000);
	// FAST - TEST
	// this->growth_rate = unifRandInt(2500,4000);
	// REAL - UNIFORM
	//this->growth_rate = unifRandInt(5000,30000);
	// REAL - WUS - BASED

	if(this->wuschel < 12){
		this->growth_rate = unifRandInt(2000,3000);
	}
	else if((this->wuschel >= 12) &&(this->wuschel <24)) {
		this->growth_rate = unifRandInt(3000,4000);
	}
	else if((this->wuschel >= 24) && (this->wuschel <36)){
		this->growth_rate = unifRandInt(4000,5000);
	}
	else if ((this->wuschel >= 36) && (this->wuschel <48)){
		this->growth_rate = unifRandInt(5000,6000);
	}
	else if ((this->wuschel >= 48) && (this->wuschel < 60)){
		this->growth_rate = unifRandInt(6000,7000);
	}	
	else if ((this->wuschel >= 60) && (this->wuschel <72)){
		this->growth_rate = unifRandInt(7000,8000);
	}
	else if ((this->wuschel >= 72) && (this->wuschel < 84)){
		this->growth_rate = unifRandInt(8000,9000);
	}
	else if ((this->wuschel >= 84) && (this->wuschel < 96)){
		this->growth_rate = unifRandInt(9000,10000);
	}
	else if((this->wuschel >= 96)&&(this->wuschel < 108)) {
		this->growth_rate = unifRandInt(10000,11000);
	}
	else if((this->wuschel >= 108)&&(this->wuschel < 120)) {
		this->growth_rate = unifRandInt(11000,12000);
	}
	else if((this->wuschel >= 120)&&(this->wuschel < 132)) {
		this->growth_rate = unifRandInt(12000,13000);
	}
	else if(this->wuschel >= 132) {
		this->growth_rate = unifRandInt(15000,16000);
	}
	if(this->cytokinin > 1200){
		this->growth_rate = unifRandInt(2000,5000);
	}


	return;
}
void Cell::set_growth_direction(Coord gd){
	this->growth_direction = gd;
	return;
}
void Cell::get_Neighbor_Cells(vector<shared_ptr<Cell>>& cells) {
	cells = this->neigh_cells;
	return;
}
void Cell::set_Left_Corner(shared_ptr<Wall_Node> new_left_corner) {
	this->left_Corner = new_left_corner;
	return;
}
void Cell::set_Wall_Count(int number_nodes) {
	this->num_wall_nodes = number_nodes;
	return;
}
double Cell::compute_membr_thresh(shared_ptr<Wall_Node> current) {
	//ran simulations changing equilibrium length to 
	//introduce a growth bias and it did not have a 
	//big effect 
	//could do more sensitivity analysis with this but 
	//for now all nodes will have same equilibrium length
	/*double l_thresh = 0;
	  double theta = 0;
	  double costheta = 0;
	  double curr_len = 0;
	  double growth_len = 0;
	  Coord curr_vec;	
	  curr_vec = current->get_Left_Neighbor()->get_Location() - current->get_Location();
	  curr_len = curr_vec.length();	
	  growth_len = 1;
	  costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
	  theta = acos( min( max(costheta,-1.0), 1.0) );
	  if(this->growth_direction == Coord(0,1)){
	  if((theta < ANGLE_FIRST_QUAD) || (theta > ANGLE_SECOND_QUAD)){
	  l_thresh = Membr_Equi_Len_Long;
	  }
	  else { 
	  l_thresh = Membr_Equi_Len_Short;
	  }
	  }*/
	double l_thresh = Membr_Equi_Len_Short;

	return l_thresh;
}

double Cell::compute_k_lin(shared_ptr<Wall_Node> current) {
	//had the idea that changing linear springs would create a
	//growth direction bias
	//but linear springs did not have the larges effect
	//for now all spring have the same linear spring constant
	/*double k_lin = 0;
	  double theta = 0;
	  double costheta = 0;
	  double curr_len = 0;
	  double growth_len = 0;
	  Coord curr_vec;	
	  curr_vec = current->get_Left_Neighbor()->get_Location() - current->get_Location();
	  curr_len = curr_vec.length();	
	  growth_len = 1;
	  costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
	  theta = acos( min( max(costheta,-1.0), 1.0) );
	  k_lin = K_LINEAR_LOOSE + K_LINEAR_STIFF*(1-pow(costheta,2));*/

	double k_lin = K_LINEAR_LOOSE;

	return k_lin;
}
double Cell::compute_k_bend(shared_ptr<Wall_Node> current) {
	//coefficient of bending spring is very important
	//nodes that are parrallel to growth direction have 
	//high bending coefficient
	//nodes that are perpendicular to growth direction have 
	//low bending coefficient
	if((growth_direction == Coord(0,1)) || (growth_direction == Coord(1,0))||(growth_direction == Coord(0,0))) {
		//fine
	}
	else{
		cout << "No growth direction assigned" << endl;
		exit(1);
	}
	double k_bend = 0;

	//This code gives preferred growth direction for the side-view model
	//by taking preferred growth direction and assigning wall nodes bending stiffness
	//by their position relative to the plane.  
	//This is only uncommented for testing purposes for tensile stress algorithm.
	if((growth_direction == Coord(0,1)) || (growth_direction == Coord(1,0))) {
		double theta = 0;
		double costheta = 0;
		double curr_len = 0;
		double growth_len = 0;
		Coord curr_vec;	
		curr_vec = current->get_Left_Neighbor()->get_Location() - current->get_Location();
		curr_len = curr_vec.length();	
		growth_len = 1;
		costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
		theta = acos( min( max(costheta,-1.0), 1.0) );
		//cout << "Theta: " << theta << endl;
		if((theta < ANGLE_FIRST_QUAD) || (theta > ANGLE_SECOND_QUAD)){
			k_bend = K_BEND_STIFF;
		} else { 
			k_bend = K_BEND_LOOSE;
		}
	} else {
		k_bend = K_BEND_UNIFORM;
	}

	//For Isotropic growth, the bending springs are baseline set to uniform stiffness.
	//This is covered in the above "Else" statement.
	//k_bend = K_BEND_UNIFORM;
	//cout << "K bend: " << k_bend << endl;
	return k_bend;
}
double Cell::compute_k_bend_div(shared_ptr<Wall_Node> current) {
	//coefficient of bending spring is very important
	//nodes that are parrallel to growth direction have 
	//high bending coefficient
	//nodes that are perpendicular to growth direction have 
	//low bending coefficient
	cout << "k bend div" << endl;
	if((growth_direction == Coord(0,1)) || (growth_direction == Coord(1,0)) || (growth_direction == Coord(0,0))) {
		//fine
	}
	else{
		cout << "No growth direction assigned" << endl;
		exit(1);
	}
	double k_bend = 0;

	//Code for computing bending stiffness for wall nodes after division.
	//Since all cells are growing isotropically, this is set to be unifrom.
	/*if((growth_direction == Coord(0,1)) || (growth_direction == Coord(1,0))){
	  double theta = 0;
	  double costheta = 0;
	  double curr_len = 0;
	  double growth_len = 0;
	  Coord curr_vec;	
	  curr_vec = current->get_Left_Neighbor()->get_Location() - current->get_Location();
	  curr_len = curr_vec.length();	
	  growth_len = 1;
	  costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
	  theta = acos( min( max(costheta,-1.0), 1.0) );
	//cout << "Theta: " << theta << endl;
	if((theta < ANGLE_FIRST_QUAD_Div) || (theta > ANGLE_SECOND_QUAD_Div)){
	k_bend = K_BEND_STIFF;
	}
	else { 
	k_bend = K_BEND_LOOSE;
	}
	} else {
	k_bend = K_BEND_UNIFORM;
	}*/

	k_bend = K_BEND_UNIFORM;
	//cout << "K bend: " << k_bend << endl;
	return k_bend;
}

void Cell::update_Wall_Angles() {
	//cout << "wall angles" << endl;
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);

#pragma omp parallel for schedule(static,1)
	for(unsigned int i=0; i< walls.size();i++) {
		//cout<< "updating" <<endl;
		walls.at(i)->update_Angle();
	}
	//cout << "Success" << endl;
	return;
}
void Cell::update_Wall_Equi_Angles() {
	//cout << "equi angles" << endl;
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
#pragma omp parallel 
	{
		double theta = 0;
		double costheta = 0;
		double curr_len = 0;
		double growth_len = 0;
		Coord curr_vec;	
		//int counter = 0;
		double new_equi_angle = 0; 
		double circle_angle  = (this->num_wall_nodes-2)*pi/(this->num_wall_nodes);
#pragma omp parallel for schedule(static,1)
		for(unsigned int i = 0; i < walls.size();i++) {	
			if(this->growth_direction != Coord(0,0)){
				curr_vec = walls.at(i)->get_Left_Neighbor()->get_Location() - walls.at(i)->get_Location();
				curr_len = curr_vec.length();	
				growth_len = this->growth_direction.length();
				costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
				theta = acos( min( max(costheta,-1.0), 1.0) );
				if((theta < ANGLE_FIRST_QUAD) || (theta > ANGLE_SECOND_QUAD)){
					new_equi_angle = pi;
				}
				else{
					//counter++;
					new_equi_angle = circle_angle;

				}
			}
			else{
				new_equi_angle = circle_angle;
			}
			new_equi_angle = circle_angle;


			//this was an idea to make the round part of the cell
			//smaller in radius but not necessary
			//if(new_equi_angle != pi) {
			//	new_equi_angle =  (counter*2-2)*pi/(counter*2);
			//}

			walls.at(i)->update_Equi_Angle(new_equi_angle);
		}
	}
	return;
}
void Cell::update_Wall_Equi_Angles_Div() {
	cout << "equi angles div" << endl;
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
#pragma omp parallel 
	{
		double theta = 0;
		double costheta = 0;
		double curr_len = 0;
		double growth_len = 0;
		Coord curr_vec;	
		//int counter = 0;
		double new_equi_angle = 0; 
		double circle_angle  = (this->num_wall_nodes-2)*pi/(this->num_wall_nodes);
#pragma omp parallel for schedule(static,1)
		for(unsigned int i = 0; i < walls.size();i++) {	
			if(this->growth_direction != Coord(0,0)){
				curr_vec = walls.at(i)->get_Left_Neighbor()->get_Location() - walls.at(i)->get_Location();
				curr_len = curr_vec.length();	
				growth_len = this->growth_direction.length();
				costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
				theta = acos( min( max(costheta,-1.0), 1.0) );
				if((theta < ANGLE_FIRST_QUAD_Div) || (theta > ANGLE_SECOND_QUAD_Div)){
					new_equi_angle = circle_angle;
				}
				else{
					//counter++;
					new_equi_angle = circle_angle;
				}
			}
			else {
				new_equi_angle = circle_angle;
			}
			//this was an idea to make the round part of the cell
			//smaller in radius but not necessary
			//if(new_equi_angle != pi) {
			//	new_equi_angle =  (counter*2-2)*pi/(counter*2);
			//}

			walls.at(i)->update_Equi_Angle(new_equi_angle);
		}
	}
	return;
}

//Calculates cell center by computing average vector of wall nodes.
void Cell::update_Cell_Center() {
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);

	Coord total_location = Coord();
#pragma omp parallel
	{
		Coord curr_loc;
#pragma omp declare reduction(+:Coord:omp_out+=omp_in) initializer(omp_priv(omp_orig))
#pragma omp for reduction(+:total_location) schedule(static,1)
		for(unsigned int i=0;i<walls.size();i++) {
			curr_loc = walls.at(i)->get_Location();
			total_location += curr_loc;
		}
	}
	this->cell_center = total_location*((1.0/static_cast<double>(num_wall_nodes)));
	//cout << cell_center<<"center"<< endl;
	return;
}

int Cell::get_Ancestry() { 
	return ancestry; 
}

void Cell::set_Div_Plane(double X, double Y) { 
	divPlaneX = X;
	divPlaneY = Y;
	return;
}

Coord Cell::get_Div_Plane(){
	Coord divPlane(divPlaneX,divPlaneY); 
	return divPlane;
}

void Cell::update_Linear_Bending_Springs(){
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	//double new_damping = this->get_Damping();
	//double k_lin = 0;
	double k_bend = 0;
	//double l_thresh = 0;

	for(unsigned int i = 0; i < walls.size();i++) {	
		//walls.at(i)->set_Damping(new_damping);
		//walls.at(i)->set_membr_len(MembrEquLen);
		//k_lin = compute_k_lin(walls.at(i));
		k_bend = compute_k_bend(walls.at(i));
		//l_thresh = compute_k_bend(walls.at(i));
		//walls.at(i)->set_K_LINEAR(k_lin);
		walls.at(i)->set_K_BEND(k_bend);
		//walls.at(i)->set_membr_len(l_thresh);
	}
	update_Wall_Equi_Angles();
	return;
}

//=============================================================
//=========================================
// Keep Track of neighbor cells and Adhesion springs
//=========================================
//=============================================================
void Cell::update_Neighbor_Cells() {
	//clear prev vector of neigh cells
	//cout << "clear" << endl;
	neigh_cells.clear();
	//grab all cells from tissue
	//cout << "cleared" << endl;
	vector<shared_ptr<Cell>> all_Cells;
	my_tissue->get_Cells(all_Cells);

	// Empty variables for holding info about other cells
	double prelim_threshold = 20;
	shared_ptr<Cell> sp_this = shared_from_this();
	//cout << "made pointer to cell" << endl;
	// iterate through all cells
#pragma omp parallel
	{
		Coord curr_Cent;
		Coord distance;
#pragma omp for schedule(static,1)
		for (unsigned int i = 0; i < all_Cells.size(); i++) {
			shared_ptr<Cell> curr = all_Cells.at(i);
			//cout << "made pointer to current neighbor" << endl;
			if (curr != sp_this) {
				curr_Cent = curr->get_Cell_Center();
				//				cout << "got center" << endl;
				// Check if cell centers are close enough together
				distance = sp_this->cell_center - curr_Cent;
				//cout << "Distance = " << distance << endl;
				if ( distance.length() < prelim_threshold ) {
#pragma omp critical
					neigh_cells.push_back(curr);
					//cout << rank << "has neighbor" << curr->get_Rank() << endl;
				}

			}
			//else you're pointing at yourself and shouldnt do anything

		}
	}

	//cout << "Cell: " << rank << " -- neighbors: " << neigh_cells.size() << endl;

	return;
}
//each cell wall node holds a vector of adhesion 
//connections and this function clears that for
//all cell wall nodes in the cell
void Cell::clear_adhesion_vectors() {
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
#pragma omp parallel
	{	

#pragma omp for schedule(static,1)	
		for(unsigned int i=0; i< walls.size();i++) {
			walls.at(i)->clear_adhesion_vec();
		}
	}
	return;
}
//for each cell wall node on current cell 
//this function searches through all the cell wall nodes on neighboring cells
//if a cell wall node on a neighboring cell is within the ADHthresh
//this function updates adhesion vector which is a private
//member variable for each cell wall node on the current cell
//this function pushes the current wall node on the neighboring cell
//onto adhesion vector
void Cell::update_adhesion_springs() {
	//get wall nodes for this cell
	vector<shared_ptr<Wall_Node>> current_cell_walls;
	this->get_Wall_Nodes_Vec(current_cell_walls);
	vector<shared_ptr<Wall_Node>> nghbr_walls_total;
	vector<shared_ptr<Wall_Node>> nghbr_walls_current;
	//int counter;
	//get all neighboring cells to this cell
	vector<shared_ptr<Cell>> neighbors;
	this->get_Neighbor_Cells(neighbors);
	for(unsigned int i = 0; i < neighbors.size(); i++) {
		neighbors.at(i)->get_Wall_Nodes_Vec(nghbr_walls_current);
		nghbr_walls_total.insert(nghbr_walls_total.end(), nghbr_walls_current.begin(), nghbr_walls_current.end());
	}	
	//	#pragma omp parallel 
	//	{
	//		#pragma omp for schedule(static,1)
	for(unsigned int i=0; i < current_cell_walls.size(); i++) {
		//counter++;
		//cout<< counter << endl;
		//cout << "Wall node" << current_cell_walls.at(i) << endl;
		current_cell_walls.at(i)->make_connection(nghbr_walls_total);
		//cout << "connection made" << endl;
	}
	//	}
	//for all cell wall nodes
	//look at adh vector, is a nodes is in the curr cell wall
	//nodes adh vector make sure the current cell wall node
	//is in that nodes adh vector
	for(unsigned int i = 0; i < current_cell_walls.size(); i++) {
		current_cell_walls.at(i)->one_to_one_check();

	}
	return;
}
//===============================================================
//============================
//  Forces, Positioning, and Shape
//============================
//===============================================================
//calls calc_forces
void Cell::calc_New_Forces(int Ti) {
	//cout << "cyts forces" << endl;	
	vector<shared_ptr<Cyt_Node>> cyts;
	this->get_Cyt_Nodes_Vec(cyts);

#pragma omp parallel for
	for (unsigned int i = 0; i < cyts.size(); i++) {
		cyts.at(i)->calc_Forces(Ti);
	}
	//cout << "cyts done" << endl;
	//calc forces on wall nodes
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	//cout<< "walls  forces" << endl;
#pragma omp parallel
	{
#pragma omp for schedule(static,1)
		for(unsigned int i=0; i < walls.size(); i++) {
			walls.at(i)->calc_Forces(Ti);
		}	
	}

	return;
}
//calls update location
//calls update cell center
//calls update wall angles
void Cell::update_Node_Locations() {
	//update cyt nodes
	vector<shared_ptr<Cyt_Node>> cyts;
	this->get_Cyt_Nodes_Vec(cyts);
#pragma omp parallel 
	{
#pragma omp for schedule(static,1)
		for (unsigned int i = 0; i < cyts.size(); i++) {
			cyts.at(i)->update_Location();
		}	
	}

	//update wall nodes
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
#pragma omp parallel 
	{	
#pragma omp for schedule(static,1)
		for(unsigned int i=0; i< walls.size();i++) {
			//cout << "update locaation" << endl;
			walls.at(i)->update_Location();
		}
	}
	//update cell_Center
	update_Cell_Center();
	//update wall_angles
	if((this->life_length == 2000)) {
		update_Wall_Equi_Angles();
	}
	update_Wall_Angles();
	//cout << "done" << endl;
	return;
}
void Cell::compute_Shape_Tensor() {
	//Initialize accumulator variable as zeros	
	//vector<vector<double>> temp(2, vector<double>(2));
	vector<vector<double>> temp;
	vector<double> tempxy;
	tempxy.push_back(0);
	tempxy.push_back(0);
	temp.push_back(tempxy);
	temp.push_back(tempxy);
	double xc,yc, xk, yk;
	/*for (int i = 0; i < 2; i++) {
	  for(int j = 0; j < 2; j++) {
	  temp[i][j] = 0;
	  }
	  }*/
	Coord cent = this->get_Cell_Center();
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	int wallcount = this->get_wall_count();
	xc = cent.get_X();
	yc = cent.get_Y();
	for (int k = 0; k < wallcount; k++) {
		xk = walls.at(k)->get_Location().get_X();
		yk = walls.at(k)->get_Location().get_Y();
		temp[0][0] += pow(xc-xk,2);
		temp[1][0] += (xc-xk) * (yc-yk);
		temp[1][1] += pow(yc-yk,2);
	}  
	//The shape tensor matrix is symmetric.
	temp[0][1] = temp[1][0];
	//Normalize
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			temp[i][j] = temp[i][j] / static_cast<double>(wallcount);
		}
	}
	shape_tensor = temp;
	return;
}

void Cell::compute_Equi_Shape_Tensor() { 
	int count = this->get_wall_count();
	double theta = 2*pi/ static_cast<double>(count);
	//The X and Y vectors represent the artifical vertices of the circular representation
	//of the cell.
	vector<double> X;
	vector<double> Y;
	for (int i = 0; i < count; i++) {
		X.push_back(EQUI_RADIUS * cos(theta*static_cast<double>(i)));
		Y.push_back(EQUI_RADIUS * sin(theta*static_cast<double>(i)));
	}
	//Initialize accumulator variable as zeros	
	//vector<vector<double>> temp(2, vector<double>(2));
	vector<vector<double>> temp;
	vector<double> tempxy;
	tempxy.push_back(0);
	tempxy.push_back(0);
	temp.push_back(tempxy);
	temp.push_back(tempxy);
	/*for (int i = 0; i < 2; i++) {
	  for(int j = 0; j < 2; j++) {
	  temp[i][j] = 0;
	  }
	  }*/
	for (int i = 0; i < count; i++) { 
		temp[0][0] += pow(X.at(i),2);
		temp[0][1] += (X.at(i)) * (Y.at(i));
		temp[1][1] += pow(Y.at(i),2);
	}
	temp[1][0] = temp[0][1];
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			temp[i][j] = temp[i][j] / static_cast<double>(count);
		}
	}
	equi_shape_tensor = temp;
	return;
}

void Cell::compute_Stress_Tensor() {
	//Initialize accumulator variable as zeros	
	this->compute_Equi_Shape_Tensor();
	//vector<vector<double>> temp(2, vector<double>(2));
	vector<vector<double>> temp;
	vector<double> tempxy;
	tempxy.push_back(0);
	tempxy.push_back(0);
	temp.push_back(tempxy);
	temp.push_back(tempxy);
	/*for (int i = 0; i < 2; i++) {
	  for(int j = 0; j < 2; j++) {
	  temp[i][j] = 0;
	  }
	  }*/
	double traceMc0 = equi_shape_tensor[0][0]+equi_shape_tensor[1][1];
	for (int i = 0; i < 2; i++) {
		for(int j = 0; j < 2; j++) {
			temp[i][j] = SHEAR_MODULUS * (shape_tensor[i][j]-equi_shape_tensor[i][j]) / traceMc0;
		}
	}
	return;
}

Coord Cell::compute_direction_of_highest_tensile_stress(){
	//Upon starting, calculate the current tensor information.
	cout << "Calculating tensors..." << endl;
	this->compute_Shape_Tensor();
	cout << "Shape Tensor Calculated." << endl;
	this->compute_Equi_Shape_Tensor();
	cout << "Equilibrium shape Tensor Calculated." << endl;
	this->compute_Stress_Tensor();
	cout << "Stress Tensor Calculated." << endl;
	vector<shared_ptr<Wall_Node>> wall_nodes;
	this->get_Wall_Nodes_Vec(wall_nodes);
	Coord max_coord;
	double max_stress = 0; 
	bool initialized = false;
	Coord curr_coord;
	Coord cent = this->get_Cell_Center();
	double curr_stress, tempX, tempY;
	for(int i = 0; i < this->get_wall_count(); i++) {
		//Vectors tested are all vectors from center to the wall nodes.
		Coord curr_coord = wall_nodes.at(i)->get_Location() - cent;
		//Normalize - tested vectors have norm length
		curr_coord = curr_coord / curr_coord.length();
		//Matrix multiply temp as S*u
		tempX = curr_coord.get_X() * stress_tensor[0][0] + curr_coord.get_Y() * stress_tensor[0][1];
		tempY = curr_coord.get_X() * stress_tensor[1][0] + curr_coord.get_Y() * stress_tensor[1][1];
		Coord temp(tempX,tempY); 
		//Compute u^T (S*u), achieved as a dot product.
		curr_stress = curr_coord.dot(temp); 
		//Accumulate the maximal vector, with some randomness as to which maximal 
		//vector is picked.
		if (curr_stress > max_stress || !initialized) {
			max_coord = curr_coord;
			initialized = true;
		} else if (curr_stress == max_stress) {
			double X = unifRand(static_cast<double>(0),static_cast<double>(1));
			if (X < 0.5) {
				max_coord = curr_coord;
			}
		}
	}
	cout << "Completed Stressvec calculation\n";
	Coord direction_vec = max_coord;
	cout << "Division plane calculated:  <" << direction_vec.get_X() << "," << direction_vec.get_Y() << ">. " << endl;
	return direction_vec;
}

Coord Cell::compute_direction_of_smallest_plane() { 
	//Upon starting, calculate the current tensor information.
	cout << "Shape Tensor Calculating...";
	this->compute_Shape_Tensor();
	cout << "Done." << endl;
	vector<shared_ptr<Wall_Node>> wall_nodes;
	this->get_Wall_Nodes_Vec(wall_nodes);
	Coord min_coord;
	bool initialized = false; 
	Coord curr_coord;
	Coord cent = this->get_Cell_Center();
	double curr_shape, tempX, tempY;
	double min_shape = 0;
	for(int i = 0; i < this->get_wall_count(); i++) {
		//Vectors tested are all vectors from center to the wall nodes.
		Coord curr_coord = wall_nodes.at(i)->get_Location() - cent;
		//Normalize - tested vectors have norm length
		curr_coord = curr_coord / curr_coord.length();
		//Matrix multiply temp as M*u
		tempX = curr_coord.get_X() * shape_tensor[0][0] + curr_coord.get_Y() * shape_tensor[0][1];
		tempY = curr_coord.get_X() * shape_tensor[1][0] + curr_coord.get_Y() * shape_tensor[1][1];
		Coord temp(tempX,tempY); 
		//Compute u^T (M*u), achieved as a dot product.
		curr_shape = curr_coord.dot(temp); 
		if (curr_shape < min_shape || !initialized) {
			min_coord = curr_coord;
			initialized = true;
		} else if (curr_shape == min_shape) {
			double X = unifRand(static_cast<double>(0),static_cast<double>(1));
			if (X < 0.5) {
				min_coord = curr_coord;
			}
		}
	}
	Coord direction_vec = min_coord;
	cout << "Division plane calculated:  <" << direction_vec.get_X() << "," << direction_vec.get_Y() << ">. " << endl;

	return direction_vec;
}

//=====================================================================
//==========================================
// Growth of Cell
//==========================================
//=====================================================================

void Cell::update_Cell_Progress(int& Ti) {
	//update life length of the current cell
	this->update_Life_Length();
	if(this->stem == 1){
		if(this->Cell_Progress < 30){
			if((Ti % growth_rate == (growth_rate -1))){
				this->add_Cyt_Node();
				this->Cell_Progress++;
			}
		}
	} else if (this->boundary == 1) {
		if(this->Cell_Progress < 30){
			if((Ti % growth_rate == (growth_rate -1))){
				this->add_Cyt_Node();
				this->Cell_Progress++;
			}
		}
	} else if(Ti<=80000 || INDEFINITE_GROWTH) {
		if((Ti % growth_rate == (growth_rate -1))){
			this->add_Cyt_Node();
			this->Cell_Progress++;
		}
	}

	return;
}
//Checks to see if division needs to happen.  If so, pass the data regarding the division back to the 
//tissue to store for later. Otherwise do nothing. DivData is of the form
//(Timestep, Dividing cell ancestry, CenterX, Center Y, DivPlaneX, DivPlaneY).
bool Cell::division_check(vector<double>& currDivData){
	vector<shared_ptr<Cell>> neighbor_cells;
	//cout <<"Before div progress" << Cell_Progress << endl;	
	if(this-> stem ==1) {
		//do nothing
	} else if(this->boundary == 1) {
		//do nothing
	} else if(this->Cell_Progress >= 30) {

		cout << "dividing cell" << this->rank <<  endl;
		//orientation of division should be 
		//fed to the division  function here
		currDivData.push_back(this->get_Ancestry());
		this->update_Cell_Center();
		Coord divCenter = get_Cell_Center();
		currDivData.push_back(divCenter.get_X());
		currDivData.push_back(divCenter.get_Y());
		shared_ptr<Cell> new_Cell = this->division();
		currDivData.push_back(divPlaneX);
		currDivData.push_back(divPlaneY);

		cout << "division success" << endl;
		//cout << "Parent cell prog" << Cell_Progress<< endl;
		//cout << "Sister cell prog" << new_Cell->get_Cell_Progress()<< endl;
		this->my_tissue->update_Num_Cells(new_Cell);
		//setting info about new cell
		//cout << "Num cells" << this->my_tissue->get_num_cells() << endl;
		new_Cell->set_Rank(this->my_tissue->get_num_cells()-1);
		//cout << "set rank" << endl;
		cout << "Parent rank: " << this->rank << endl;
		cout << "Sister rank: " << new_Cell->get_Rank() << endl;
		cout << new_Cell->get_wall_count() << endl;
		cout << new_Cell->get_cyt_count() << endl;
		cout << this->get_wall_count() << endl;
		cout << this->get_cyt_count() << endl;
		cout << "Parent" << this << endl;
		cout << "Parent progress: " << this->get_Cell_Progress() << endl;
		cout << "New cell" << new_Cell << endl;
		cout << "New progress: " << new_Cell->get_Cell_Progress() << endl;

		//layer in division function		
		//damping in division function
		//boundary needs to be figured out
		//life length set to 0 in constructor
		//life length of parent cell reset in 
		//division function
		//cyt nodes in divison function
		//wall nodes in division function
		//all cell progress set to 0 in division 
		//cell center in division function
		//cyt and wus in division function
		//growth rate in div function
		//growth direction inherited in div function
		//left corner in divison function  
		//cout << "adhesion division" << endl;
		new_Cell->update_Neighbor_Cells();
		new_Cell->update_adhesion_springs();
		new_Cell->get_Neighbor_Cells(neighbor_cells);
		for(unsigned int i =0; i < neighbor_cells.size(); i++) {
			neighbor_cells.at(i)->update_Neighbor_Cells();
			neighbor_cells.at(i)->clear_adhesion_vectors();
			neighbor_cells.at(i)->update_adhesion_springs();

		}
		return true;
	}
	return false;
}

double Cell::calc_Area() {
	shared_ptr<Wall_Node> curr = left_Corner;
	shared_ptr<Wall_Node> next = NULL;
	shared_ptr<Wall_Node> orig = curr;
	Coord a_i;
	Coord a_j;
	double area = 0;
	double curr_area = 0;
	do {
		next = curr->get_Left_Neighbor();
		a_i = curr->get_Location() - cell_center;
		a_j = next->get_Location() - cell_center;
		curr_area = 0.5*sqrt(pow(a_i.cross(a_j),2));
		area += curr_area;
		curr = next;
	} while(next != orig);
	//cout << "Area: " << area << endl;
	return area;
}
void Cell::add_wall_Node_Check(int Ti) {
	//cout << "adding a wall node" << endl;
	add_Wall_Node(Ti);
	return;
}
void Cell::delete_wall_Node_Check(int Ti){
	delete_Wall_Node(Ti);
	return;
}
void Cell::add_Wall_Node(int Ti) {

	//find node to the right of largest spring
	shared_ptr<Cell> this_cell= shared_from_this();
	shared_ptr<Wall_Node> right;
	//cout  << "Find largest length" << endl;
	find_Largest_Length(right);
	//cout << "Largest found" << endl;
	shared_ptr<Wall_Node> left;
	Coord location;
	double l_thresh;
	double k_lin;
	double k_bend;
	if(this->life_length>2000){	
		if(right != NULL) {
			//find location and set neighbors for new node
			left = right->get_Left_Neighbor();
			location  = (right->get_Location() + left->get_Location())*0.5;
			shared_ptr<Wall_Node> added_node = make_shared<Wall_Node>(location, this_cell, left, right);
			this->add_wall_node_vec(added_node);
			double new_damping = this->get_Damping();
			right->set_Left_Neighbor(added_node);
			left->set_Right_Neighbor(added_node);
			//set the variables for the new node
			l_thresh = compute_membr_thresh(added_node);
			k_lin = compute_k_lin(added_node);
			k_bend = compute_k_bend(added_node);
			added_node->set_Damping(new_damping);
			added_node->set_K_LINEAR(k_lin);
			added_node->set_K_BEND(k_bend);
			added_node->set_membr_len(l_thresh);	
			added_node->set_added(1);
			//adhesion for the new node
			vector<shared_ptr<Cell>> neighbors;
			this->get_Neighbor_Cells(neighbors);
			//cout << "adh added node find closest" << endl;
			vector<shared_ptr<Wall_Node>> nghbr_walls_total;
			vector<shared_ptr<Wall_Node>> nghbr_walls_current;
			for(unsigned int i = 0; i < neighbors.size(); i++) {
				neighbors.at(i)->get_Wall_Nodes_Vec(nghbr_walls_current);
				nghbr_walls_total.insert(nghbr_walls_total.end(), nghbr_walls_current.begin(), nghbr_walls_current.end());
			}
			added_node->make_connection(nghbr_walls_total);
			//cout << "adh added node success" << endl;
			//update angles
			//should i update k_bennnnndddddd?????
			update_Wall_Equi_Angles();
			update_Wall_Angles();
		}
	}
	return;
}
void Cell::delete_Wall_Node(int Ti) {
	shared_ptr<Wall_Node> left = NULL;
	shared_ptr<Wall_Node> right = NULL;
	shared_ptr<Wall_Node> small = NULL;
	//vector<Cell*>neighbors;

	this->find_Smallest_Length(small);
	if(small !=NULL) {
		cout << "delete initiated" << endl;
		left = small->get_Left_Neighbor();
		right = small->get_Right_Neighbor();

		//if small is the left corner cell reassign
		if(this->left_Corner == small) {
			//cout << " set left corner" << endl;
			this->set_Left_Corner(left);
		}
		//need to make sure all nodes connected to small
		//via adhesion are erased
		//small->remove_from_adh_vecs();
		//small->clear_adh_vec();

		//set new neighbors so nothing points at small
		left->set_Right_Neighbor(right);
		right->set_Left_Neighbor(left);
		//reindex
		this->wall_nodes.clear();
		this->num_wall_nodes = 0;

		shared_ptr<Wall_Node> curr = this->left_Corner;
		shared_ptr<Wall_Node> next = NULL;
		shared_ptr<Wall_Node> orig = curr;

		do {
			this->wall_nodes.push_back(curr);
			next = curr->get_Left_Neighbor();
			num_wall_nodes++;
			curr = next;
		}while(next != orig);
		//cout << "update equi angles" << endl;

		update_Wall_Equi_Angles();

		//cout << "update angles" << endl;
		update_Wall_Angles();
	}
	return;
}
//finds right neighbor node of smallest length on membrane
void Cell::find_Smallest_Length(shared_ptr<Wall_Node>& right) {
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	double max_len = 100;
	//#pragma omp parallel
	//{
	shared_ptr<Wall_Node> left_neighbor;
	double curr_len = 0;
	//#pragma omp for schedule(static,1)
	for (unsigned int i = 0; i < walls.size();i++) {
		left_neighbor = walls.at(i)->get_Left_Neighbor();
		curr_len = (walls.at(i)->get_Location()-left_neighbor->get_Location()).length();
		if(curr_len < .05){
			if(curr_len < max_len) {
				//#pragma omp critical
				max_len = curr_len;
				right = walls.at(i);
			}
		}
	}
	//}
	return;
}
//finds right neighbor node of largest length on membrane
void Cell::find_Largest_Length(shared_ptr<Wall_Node>& largest) {
	vector<shared_ptr<Wall_Node>> walls; 
	this->get_Wall_Nodes_Vec(walls);
	double max_len = 0;
	shared_ptr<Wall_Node> right;
	shared_ptr<Wall_Node> current;
	// Unused double theta = 0;
	// Unused double costheta = 0;
	// Unused double curr_len = 0;
	// Unused double growth_len = 0;
	Coord curr_vec;	

	//#pragma omp parallel 
	shared_ptr<Wall_Node> left_neighbor;
	double current_len = 0;
	//#pragma omp for schedule(static,1)
	for(unsigned int i=0; i < walls.size(); i++) {
		current = walls.at(i);
		curr_vec = current->get_Left_Neighbor()->get_Location() - current->get_Location();
		//curr_len = curr_vec.length();	
		// Unused growth_len = 1;
		// Unused costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
		// Unused theta = acos( min( max(costheta,-1.0), 1.0) );

		//if((theta < ADD_WALL_NODE_ANGLE_FIRST_QUAD) ||(theta > ADD_WALL_NODE_ANGLE_SECOND_QUAD))
		left_neighbor = walls.at(i)->get_Left_Neighbor();
		current_len = (walls.at(i)->get_Location()-left_neighbor->get_Location()).length();
		//if(curr_len > MEMBR_THRESH_LENGTH) 
		if(current_len > max_len){
			max_len = current_len;
			right = walls.at(i);


		}	

	}
	largest = right;

	return;

}

void Cell::add_Cyt_Node() {
	double new_damping = this->get_Damping();
	shared_ptr<Cell> this_cell = shared_from_this();
	shared_ptr<Cyt_Node> cyt = make_shared<Cyt_Node>(cell_center, this_cell);
	cyt_nodes.push_back(cyt);
	cyt->set_Damping(new_damping);

	num_cyt_nodes++;
	return;
}
//===========================================================
//==================================
// Output Functions
//==================================
//===========================================================
void Cell::print_Data_Output(ofstream& ofs) {
	ofs << "This is where data output goes" << endl;

	return;
}

int Cell::update_VTK_Indices(int& id) {
	//cout << "ID before: " << id << endl;
	int rel_cnt = 0;

	shared_ptr<Wall_Node> curr_wall = left_Corner;
	do { 
		curr_wall->update_VTK_Id(id);
		id++;
		for(unsigned int i = 0; i < curr_wall->get_adh_vec().size(); i++){
			rel_cnt++;
		}
		curr_wall = curr_wall->get_Left_Neighbor();
	} while (curr_wall != left_Corner);

	for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
		cyt_nodes.at(i)->update_VTK_Id(id);
		id++;
	}
	//cout << "ID after: " << id << endl;
	return rel_cnt;
}
void Cell::print_VTK_Adh(ofstream& ofs) {
	int my_id, nei_id;
	shared_ptr<Wall_Node> neighbor = NULL;
	shared_ptr<Wall_Node> curr_wall = left_Corner;
	vector<shared_ptr<Wall_Node>> nodes;
	do {
		for(unsigned int i = 0; i < curr_wall->get_adh_vec().size(); i++){
			nodes = curr_wall->get_adh_vec();
			neighbor = nodes.at(i);
			if(neighbor != NULL) {
				my_id = curr_wall->get_VTK_Id();
				nei_id = neighbor->get_VTK_Id();
				ofs << 2 << ' ' << my_id << ' ' << nei_id << endl;
			}
		}
		curr_wall = curr_wall->get_Left_Neighbor();
	} while(curr_wall != left_Corner);
	return;
}
Coord Cell::average_coordinates(){
	Coord direction = Coord(0,0);
	for(unsigned int i = 0; i <wall_nodes.size() ;i++){
		direction = direction + wall_nodes.at(i)->get_Location();
	}
	double num_nodes = (double) wall_nodes.size();
	direction = direction/num_nodes;
	return direction;
}
void Cell::print_direction_vec(ofstream& ofs){
	//average coordinates of all nodes
	Coord sum = this->average_coordinates();
	//point 1 is center + 1*direction vector
	Coord point1 = this->cell_center + sum*.05;
	Coord point2 = this->cell_center - sum*.05;
	ofs << point1.get_X() << ' ' << point1.get_Y() << ' ' << 0 << endl;
	ofs << point2.get_X() << ' ' << point2.get_Y() << ' ' << 0 << endl;

	return;
}
void Cell::print_locations(ofstream& ofs) {
	ofs << this->get_Rank() << endl;
	shared_ptr<Wall_Node> curr_wall = left_Corner;
	shared_ptr<Wall_Node> orig = curr_wall;
	//	cout << "knows left corner" << endl;
	do {


		Coord loc = curr_wall->get_Location();
		ofs << loc.get_X() << ' ' << loc.get_Y() << ' ' << 0 <<' '<< 1 << endl;
		//cout<< "maybe cant do left neighbor" << endl;
		curr_wall = curr_wall->get_Left_Neighbor();

		//cout << "did it  " << count << endl;
	} while (curr_wall != orig);

	//cout << "walls worked" << endl;
	/*for (unsigned int i = 0; i < cyt_nodes.size(); i++) {

	  Coord loc = cyt_nodes.at(i)->get_Location();
	  ofs << loc.get_X() << ' ' << loc.get_Y() << ' ' << 0 << 0 << endl;
	  }*/
	return;
}

void Cell::print_VTK_Points(ofstream& ofs, int& count) {
	shared_ptr<Wall_Node> curr_wall = left_Corner;
	shared_ptr<Wall_Node> orig = curr_wall;
	//	cout << "knows left corner" << endl;
	do {
		Coord loc = curr_wall->get_Location();
		ofs << loc.get_X() << ' ' << loc.get_Y() << ' ' << 0 << endl;
		//cout<< "maybe cant do left neighbor" << endl;
		curr_wall = curr_wall->get_Left_Neighbor();
		count++;
		//cout << "did it  " << count << endl;
	} while (curr_wall != orig);

	//cout << "walls worked" << endl;
	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		Coord loc = cyt_nodes.at(i)->get_Location();
		ofs << loc.get_X() << ' ' << loc.get_Y() << ' ' << 0 << endl;
		count++;
	};

	//cout << "points worked" << endl;
	return;
}

void Cell::print_VTK_Scalars_Average_Pressure(ofstream& ofs) {
	//float pressure = this->average_Pressure();
	shared_ptr<Wall_Node> curr_wall = left_Corner;
	do {
		//concentration = curr_wall->get_My_Cell()->get_WUS_concentration();
		//	ofs << pressure << endl;

		curr_wall = curr_wall->get_Left_Neighbor();
	} while (curr_wall != left_Corner);
	//for(unsigned int i=0;i<wall_nodes.size();i++){
	//	ofs << pressure << endl;
	//}
	for(unsigned int i=0;i<cyt_nodes.size();i++){
		//	ofs<< pressure << endl;
	}
	return;
}

void Cell::print_VTK_Scalars_WUS(ofstream& ofs) {

	double concentration = 0;
	shared_ptr<Wall_Node> curr_wall = left_Corner;
	do {
		concentration = curr_wall->get_My_Cell()->get_WUS_concentration();
		ofs << concentration << endl;

		curr_wall = curr_wall->get_Left_Neighbor();

	} while (curr_wall != left_Corner);


	for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
		concentration = cyt_nodes.at(i)->get_My_Cell()->get_WUS_concentration();
		ofs << concentration << endl;
	}
	return;
}
void Cell::print_VTK_Scalars_CK(ofstream& ofs) {

	double concentration = 0;
	shared_ptr<Wall_Node> curr_wall = left_Corner;
	do {
		concentration = curr_wall->get_My_Cell()->get_CYT_concentration();
		ofs << concentration << endl;

		curr_wall = curr_wall->get_Left_Neighbor();

	} while (curr_wall != left_Corner);


	for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
		concentration = cyt_nodes.at(i)->get_My_Cell()->get_CYT_concentration();
		ofs << concentration << endl;
	}
	return;
}
void Cell::print_VTK_Scalars_Node(ofstream& ofs) {
	shared_ptr<Wall_Node> currW = left_Corner;
	double color;
	do {
		if(currW->get_added()==1){
			color = 30.0;
		}
		else{
			color = 0.0;
		}
		ofs << color << endl;
		currW = currW->get_Left_Neighbor();
	} while(currW != left_Corner);
	for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
		color = 0.0;
		ofs << color << endl;
	}
	return;
}



























///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////OLD CODE////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

/*Coord Cell::compute_direction_of_highest_tensile_stress(){
//average position of all cell wall nodes
vector<shared_ptr<Wall_Node>> wall_nodes;
this->get_Wall_Nodes_Vec(wall_nodes);
Coord next_coord;
Coord curr_coord;
Coord direction_vec(0,0);
shared_ptr<Wall_Node> curr = wall_nodes.at(0);
shared_ptr<Wall_Node> next = curr;
//double min_dist = (curr->get_Location() - (curr->get_Left_Neighbor())->get_Location()).length();
//double curr_dist;
int counter = 0;
cout << "Tensile_stress Do-while start" << endl;
//Distance is equivalent to stress due to constant equilibrium dist.
do{
curr = next;
next = curr->get_Left_Neighbor();
curr_coord = curr->get_Location();
next_coord = next->get_Location();
//curr_dist = (curr_coord - next_coord).length();
//min_dist = ((min_dist >= curr_dist) ? curr_dist : min_dist);
counter++;
} while(next != wall_nodes.at(0));
cout << "Tensile_stress Do-while end, Counter = " << counter << endl;
// Minimium distance between nodes should never be 0.
//CHECK
counter = this->get_wall_count();*/
/*if (min_dist == 0) { 
  std::cout << "min_dist = 0, nodes are overlapped" << endl;
  exit(1);
  } */ 
/*//Initialize arrays for each node's stress direction and its relative tensile stress
  Coord * stress_vec =  new Coord [counter];
  double * stress_mag = new double [counter];
  for (int i = 0; i < counter; i++) {
  stress_mag[i] = 0;
  }
  shared_ptr<Wall_Node> temp = (curr->get_Left_Neighbor())->get_Left_Neighbor();
//Calculate relative stresses based on distance stretched. 
for (int j = 0; j < counter; j++) {
for (int i = 0; i < 4; i++) {
stress_vec[j] = (curr->get_Location() - this->get_Cell_Center()).perpVector();
//stress_mag[j] = stress_mag[j] + ((temp->get_Right_Neighbor())->get_Location() - temp->get_Location()).length() - min_dist;
stress_mag[j] = stress_mag[j] + ((temp->get_Right_Neighbor())->get_Location() - temp->get_Location()).length() - Membr_Equi_Len_Long;
temp = temp->get_Right_Neighbor();
}
}
double sumstress = 0; 
for (int i = 0; i < counter; i++) {
sumstress = sumstress + stress_mag[i];
}
//Calculate Average stresses and average out the fectors.
int alignment;
//Calculates the weighted average of vectors.  Alignment requires that 
//the vectors are <= 90 degrees away from one another, otherwise we may find
//vectors x,y with angle as ~179 degrees as cancelling out instead of contributing
//toward the same direction of growth.
for (int i = 0; i < counter; i++) {
alignment = (direction_vec.dot(stress_vec[i]) >= 0) ? 1 : -1; 
direction_vec = direction_vec + stress_vec[i]*alignment*(stress_mag[i] / sumstress);
if (direction_vec.get_Y() < 0) {
direction_vec = direction_vec * (-1); 
}
}

delete stress_vec;
delete stress_mag; 
cout << "Completed Stressvec calculation\n";
return direction_vec;
}*/

/*void Cell::print_VTK_Scalars_Wall_Pressure(ofstream& ofs){
  this->wall_Pressure();
  Wall_Node* curr_wall = this->left_Corner;

  double pressure;
  do{
  pressure = curr_wall->get_pressure();
  curr_wall = curr_wall->get_Left_Neighbor();
  ofs << pressure << endl;
  }while(curr_wall!= left_Corner);
  for(unsigned int i=0; i<cyt_nodes.size(); i++){
  ofs << 0 << endl;
  }
  return;
  }*/
/*void Cell::print_VTK_Scalars_Average_Pressure_cell(ofstream& ofs) {
  double pressure = this->average_Pressure();
  for(unsigned int i=0;i<1;i++){
  ofs << pressure << endl;
  }
  return;
  }*/

/*void Cell::print_VTK_Scalars_WUS_cell(ofstream& ofs) {
  for(unsigned int i = 0; i < 1; i++) {

  double concentration = cyt_nodes.at(i)->get_My_Cell()->get_WUS_concentration();

  ofs << concentration << endl;
  }
  return;
  }

  void Cell::print_VTK_Scalars_CYT(ofstream& ofs) {

//	double concentration = 0;

shared_ptr<Wall_Node> curr_wall = left_Corner;
do {
double concentration = curr_wall->get_My_Cell()->get_CYT_concentration();
ofs << concentration << endl;

curr_wall = curr_wall->get_Left_Neighbor();

} while (curr_wall != left_Corner);


for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
double concentration = cyt_nodes.at(i)->get_My_Cell()->get_CYT_concentration();
ofs << concentration << endl;
}
return;
}*/
/*Coord Cell::compute_point_on_line(double t){
  Coord r_0 = this->get_Cell_Center();
  Coord v = this->compute_longest_axis();
  Coord p = Coord(r_0.get_X() + v.get_X()*t, r_0.get_Y() + v.get_Y()*t);
  return p;
  }*/


/*Coord Cell::compute_direction_of_highest_tensile_stress(){
//average position of all cell wall nodes
vector<shared_ptr<Wall_Node>> wall_nodes;
this->get_Wall_Nodes_Vec(wall_nodes);
Coord next_coord;
Coord curr_coord;
Coord direction_vec;
shared_ptr<Wall_Node> curr = wall_nodes.at(0);
shared_ptr<Wall_Node> next;
double delta_x = 0;
double delta_y = 0;
double x = 0;
double y = 0;
double average_x;
double average_y;
int counter = 0;
//Distance is equivalent to stress due to constant equilibrium dist.
do{
next = curr->get_Left_Neighbor();
curr_coord = curr->get_Location();
next_coord = next->get_Location();
delta_x = next_coord.get_X() - curr_coord.get_X();
delta_y = next_coord.get_Y() - curr_coord.get_Y();
x = x + delta_x;
y = y + delta_x;
counter++;
} while(next != wall_nodes.at(0));
average_x = x/counter;
average_y = y/counter;

direction_vec = Coord(-average_y,average_x);
return direction_vec;
}*/

/*		   
		   void Cell::print_VTK_Scalars_Total(ofstream& ofs) {

//	double concentration = 0;

Wall_Node* curr_wall = left_Corner;
do {
double concentration = curr_wall->get_My_Cell()->get_total_concentration();
ofs << concentration << endl;

curr_wall = curr_wall->get_Left_Neighbor();

} while (curr_wall != left_Corner);


for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
double concentration = cyt_nodes.at(i)->get_My_Cell()->get_total_concentration();
ofs << concentration << endl;
}
return;
}

void Cell::print_VTK_Vectors(ofstream& ofs) {

Wall_Node* curr_wall = left_Corner;
do {
Coord force = curr_wall->get_CytForce();
ofs << force.get_X() << ' ' << force.get_Y() << ' ' << 0 << endl;

curr_wall = curr_wall->get_Left_Neighbor();

} while(curr_wall != left_Corner);

for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
Coord force = cyt_nodes.at(i)->get_Force();
}
return;
}

void Cell::print_VTK_Vectors(ofstream& ofs) {

Wall_Node* curr_wall = left_Corner;
do {
Coord force = curr_wall->get_CytForce();
ofs << force.get_X() << ' ' << force.get_Y() << ' ' << 0 << endl;

curr_wall = curr_wall->get_Left_Neighbor();

} while(curr_wall != left_Corner);

for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
Coord force = cyt_nodes.at(i)->get_Force();
ofs << force.get_X() << ' ' << force.get_Y() << ' ' << 0 << endl;
}

return;
}*/





//////////////////////////////////
