
// Includes
#include <vector>
#include <cmath>
#include <ctime>
#include <exception>
#include <iostream>

#include "ros/ros.h"
#include "geometry_msgs/Twist.h"
#include "geometry_msgs/PoseArray.h"
#include "nav_msgs/Odometry.h"
#include "tf/transform_datatypes.h"
#include "eigen3/Eigen/Dense"
#include "/home/johan/catkin_ws/src/control_project/control_controller/include/eiquadprog.hpp"

#define _USE_MATH_DEFINES

namespace mpc {

  struct Tools {
    static Eigen::Matrix3d computeStateMatrix(
      const Eigen::Vector3d& x, 
      double wheel_radius, 
      double old_linear_vel)
    {
      Eigen::Matrix3d A;
      A <<  0.0, 0.0, (-wheel_radius * old_linear_vel * std::sin(x[2])),
            0.0, 0.0, (wheel_radius * old_linear_vel * std::cos(x[2])),
            0.0, 0.0, 0.0;
      return A;
    }

    static Eigen::Matrix<double, 3, 2> computeControlMatrix(
      const Eigen::Vector3d& x,
      double wheel_radius,
      double base_lenght)
    {
      Eigen::Matrix<double, 3, 2> B;
      B <<  (wheel_radius * std::cos(x[2])), 0.0,
            (wheel_radius * std::sin(x[2])), 0.0,
            0.0, (wheel_radius/base_lenght);
      return B;
    }

    static Eigen::Vector3d computeXdot(
      const Eigen::Vector3d& x, 
      const Eigen::Vector3d& next_x, 
      double dt)
    {
      Eigen::Vector3d x_dot;
      x_dot <<  ((next_x[0] - x[0]) / dt),
                ((next_x[1] - x[1]) / dt),
                ((next_x[2] - x[2]) / dt);
      return x_dot;
    }

    static double subtractAngles(double angle1, double angle2){
      double angle = angle1 - angle2;
      
      if(angle > M_PI)
        return angle - 2.0 * M_PI;
      else if(angle < -M_PI)
        return angle + 2.0 * M_PI;
      return angle;
    }
  };

  class ControlVectors {

    std::vector<Eigen::Vector2d> ref_u;
    double _wheel_radius;
    double _base_lenght;
    double _old_linear_vel;

    public:

    ControlVectors(double wheel_radius = 1.65, double base_lenght = 2.5){
      this->_wheel_radius = wheel_radius;
      this->_base_lenght = base_lenght;
      this->_old_linear_vel = 0.0;
    }

    void push_back(const Eigen::Vector3d& x, const Eigen::Vector3d& next_x, double dt){
      
      Eigen::Matrix3d A = Tools::computeStateMatrix(x, _wheel_radius, _old_linear_vel);
      Eigen::Matrix<double, 3, 2> B = Tools::computeControlMatrix(x, _wheel_radius, _base_lenght);
      Eigen::Vector3d x_dot = Tools::computeXdot(x, next_x, dt);
      
      Eigen::Vector3d sum = x_dot - A * x;
      Eigen::Vector2d u = B.colPivHouseholderQr().solve(sum);

      this->_old_linear_vel = u[0];
      this->ref_u.push_back(u);
    }

    void push_back(double linear_vel, double angular_vel){
      Eigen::Vector2d u;
      u << linear_vel, angular_vel;
      this->ref_u.push_back(u);
    }

    void push_back(Eigen::Vector2d vector){
      this->ref_u.push_back(vector);
    }

    void pop_back(){ ref_u.pop_back(); }

    void clear(){ ref_u.clear(); }

    Eigen::Vector2d &operator[](int n){ return ref_u[n]; }
    
    const Eigen::Vector2d &operator[](int n) const { return ref_u[n]; }

    Eigen::Vector2d get(int n) const { return ref_u[n]; }

    std::size_t size() const { return ref_u.size(); }

    ControlVectors* getHorizon(int start, int length){
      if(start+length+1 > ref_u.size() || start < 0)
        throw std::out_of_range::exception();

      ControlVectors* horizon = new ControlVectors();
      for(int i = start; i < (start + length); i++){
        horizon->push_back(ref_u[i]);
      }
      return horizon;
    }

    ~ControlVectors(){}
  };

  class StateVectors {

    std::vector<Eigen::Vector3d> ref_x;

    public:

    StateVectors() {} 

    void push_back(const geometry_msgs::Pose& pose){
      double x = pose.position.x;
      double y = pose.position.y;
      double theta = tf::getYaw(pose.orientation);
      this->ref_x.push_back(Eigen::Vector3d(x, y, theta));
    }

    void push_back(Eigen::Vector3d vector){
      this->ref_x.push_back(vector);
    }

    void pop_back(){ ref_x.pop_back(); }

    void clear(){ ref_x.clear(); }

    std::size_t size() const { return ref_x.size(); } 

    Eigen::Vector3d &operator[](int n){ return ref_x[n]; }

    const Eigen::Vector3d &operator[](int n) const { return ref_x[n]; }

    Eigen::Vector3d get(int n) const { return ref_x[n]; }

    StateVectors* getHorizon(int start, int length){
      if(start+length+1 > ref_x.size() || start < 0)
        throw std::out_of_range::exception();

      StateVectors* horizon = new StateVectors();
      for(int i = start; i < (start + length); i++){
        horizon->push_back(ref_x[i]);
      }
      return horizon;
    }

    ~StateVectors(){}
  };

  class MPC {

    private:

    double _q, _r;
    double _wheel_radius;
    double _base_lenght;
    int _horizon_lenght;

    Eigen::Matrix<double,-1,3> _S;
    Eigen::MatrixXd _T;
    Eigen::MatrixXd _P;
    Eigen::VectorXd _p;

    public:

    MPC(double q = 10.0, double r = 0.1, double wheel_radius = 1.0,  double base_length = 2.5){
      this->_q = q;
      this->_r = r;
      this->_wheel_radius = wheel_radius;
      this->_base_lenght = base_length;
    }

    void setStateWeight(double q){ this->_q = q; }

    void setControlWeight(double r){ this->_r = r; }

    void setWheelRadius(double value){ this->_wheel_radius = value; }

    void setBaseLength(double value){ this->_base_lenght = value; }

    double getStateWeight(){ return this->_q; }

    double getControlWeight(){ return this->_r; }

    double getWheelRadius(){ return this->_wheel_radius; }

    double getBaseLength(){ return this->_base_lenght; }

    Eigen::Vector2d computeControl(
      const StateVectors* horizon_x_ref, 
      const ControlVectors* horizon_u_ref, 
      Eigen::Vector3d current_x, 
      double dt)
    {
      if(horizon_x_ref->size() != horizon_u_ref->size()){
        throw std::length_error::exception();
      }

      _horizon_lenght = (int)horizon_x_ref->size();
      ROS_INFO("Horizon length = %d", _horizon_lenght);
      
      computeST(horizon_x_ref, horizon_u_ref, dt);
      ROS_INFO("S=[%d,%d], T=[%d,%d]", (int)_S.rows(), (int)_S.cols(), (int)_T.rows(), (int)_T.cols());

      Eigen::Vector3d x_zero = current_x - horizon_x_ref->get(0);

      int D_R = _horizon_lenght * 2;
      int D_Q = _horizon_lenght * 3;
      Eigen::MatrixXd R = Eigen::MatrixXd::Identity(D_R, D_R) * _r;
      Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(D_Q, D_Q) * _q;
      ROS_INFO("R=[%d,%d], Q=[%d,%d]", (int)R.rows(), (int)R.cols(), (int)Q.rows(), (int)Q.cols());

      computePp(R, Q, x_zero);
      ROS_INFO("P=[%d,%d], p=[%d,%d]", (int)_P.rows(), (int)_P.cols(), (int)_p.rows(), (int)_p.cols());

      // Quadratic programming and constraints
      Eigen::MatrixXd CE = Eigen::MatrixXd::Zero(_P.rows(), _P.cols());
      Eigen::VectorXd ce0 = Eigen::VectorXd::Zero(_P.rows());

      Eigen::MatrixXd CI = Eigen::MatrixXd::Zero(_P.rows(), _P.cols());
      Eigen::VectorXd ci0 = Eigen::VectorXd::Zero(_P.rows());

      Eigen::MatrixXd double_P = _P * 2.0;
      Eigen::VectorXd u_min; 
      Eigen::solve_quadprog(double_P, _p, CE, ce0, CI, ci0, u_min);
      ROS_INFO("u_vec=[%d,%d]", (int)u_min.rows(), (int)u_min.cols());

      Eigen::Vector2d u;
      u << u_min[0], u_min[1]; 
      return u + horizon_u_ref->get(0);
    }

    ~MPC(){}

    private:
    
    void computeST(
      const StateVectors* horizon_x_ref, 
      const ControlVectors* horizon_u_ref,
      double dt)
    {
      _S = Eigen::Matrix<double,-1,3>::Zero(_horizon_lenght*3, 3);
      _T = Eigen::MatrixXd::Zero(_horizon_lenght*3, _horizon_lenght*2);

      Eigen::Matrix3d SA;
      double old_linear_vel = 0.0;
      for(int i = 0; i < _horizon_lenght; i++){
        Eigen::Matrix3d A = Tools::computeStateMatrix(horizon_x_ref->get(i), _wheel_radius, old_linear_vel) * dt + Eigen::Matrix3d::Identity();
        
        if(i > 0)
          SA = SA * A;
        else
          SA = A;

        _S.block(i*3, 0, 3, 3) = SA;
       
        Eigen::Matrix<double, 3, 2> B = Tools::computeControlMatrix(horizon_x_ref->get(i), _wheel_radius, _base_lenght) * dt;
        _T.block(i*3, i*2, 3, 2) = B;
        if(i > 0){
          for(int j = 0; j < i; j++){
            _T.block(i*3, j*2, 3, 2) = A * _T.block((i-1)*3, j*2, 3, 2);
          }
        }
        
        old_linear_vel = (horizon_u_ref->get(i))[0];
      }
    }

    void computePp(
      Eigen::MatrixXd& R, 
      Eigen::MatrixXd& Q, 
      Eigen::Vector3d& x_zero)
    {
      Eigen::MatrixXd T_transpose = _T;
      T_transpose.transposeInPlace();
      
      _P = T_transpose * Q * _T + R;
      _p = 2*T_transpose * Q * _S * x_zero;
    }
  };

};

class Pose2D{
    
    public:

    Pose2D(double x = 0.0, double y = 0.0, double theta = 0.0){
      this->_pose << x, y, theta; 
    }

    double getX() { return _pose[0]; }

    void setX(double x) { this->_pose[0] = x; }

    double getY() { return _pose[1]; }

    void setY(double y) { this->_pose[1] = y; }

    double getTheta() { return _pose[2]; }

    void setTheta(double theta) { this->_pose[2] = theta; }

    Eigen::Vector3d getPoseVector() { return _pose; } 

    private:

    Eigen::Vector3d _pose;
};

class Robot {

  ros::NodeHandle                   nh;
  ros::Subscriber                   traj_sub;
  ros::Subscriber                   odom_sub;
  ros::Publisher                    movement_pub;

  std::vector<geometry_msgs::Pose>  trajectory;
  mpc::StateVectors*                ref_state;
  mpc::ControlVectors*              ref_control;
  mpc::MPC*                         controller;
  Pose2D                            pose;
  
  enum RunState{
    Init,
    Running,
    Idle
  };
  
  RunState                          state;
  Eigen::Vector2d                   cmd_u;
  bool                              traj_recived;
  bool                              run_init;
  int                               run_counter;
  struct timespec                   time_now;
  struct timespec                   time_since_last;
  double                            dt;
  unsigned long                     nano_dt;

  void trajectoryCallback(const geometry_msgs::PoseArray::ConstPtr& msg) {
    this->trajectory.clear();
    this->trajectory = msg->poses;
    ROS_INFO("Trajectory size %d", (int)trajectory.size());

    ref_state->clear();
    for(auto traj : this->trajectory){
      ref_state->push_back(traj);
    }
    ROS_INFO("Reference state size %d", (int)ref_state->size());

    ref_control->clear();
    for(int i = 0; i < (ref_state->size()-1); i++){
      ref_control->push_back(ref_state->get(i), ref_state->get(i+1), dt);
    }
    ref_control->push_back(0.0, 0.0);
    ROS_INFO("Control state size %d", (int)ref_control->size());

    state = RunState::Init;
  }

  void odomCallback(const nav_msgs::Odometry::ConstPtr& msg) {

    geometry_msgs::Pose robot_pose = msg->pose.pose;
    this->pose.setX(robot_pose.position.x);
    this->pose.setY(robot_pose.position.y);
    this->pose.setTheta(tf::getYaw(robot_pose.orientation));
    //ROS_INFO("Recived new pose [%.2f, %.2f, %.2f]", pose.getX(), pose.getY(), pose.getTheta());
  }

  public:

  Robot(double wheel_radius, double base_length, double dt) {
    this->movement_pub = nh.advertise<geometry_msgs::Twist>("/robot_0/cmd_vel", 10);    
    this->traj_sub = nh.subscribe("/planner/trajectory", 10, &Robot::trajectoryCallback, this);
    this->odom_sub = nh.subscribe("/robot_0/base_pose_ground_truth", 100, &Robot::odomCallback, this);
    
    this->ref_state = new mpc::StateVectors();
    this->ref_control = new mpc::ControlVectors(wheel_radius, base_length);
    this->controller = new mpc::MPC(1.0, 0.1, wheel_radius, base_length);

    //Seconds to nanoseconds
    this->dt = dt;
    this->nano_dt = (int)(dt * 1000000000.0);
    this->cmd_u = Eigen::Vector2d::Zero();
    this->state = RunState::Idle;
  }

  ~Robot(){
    delete ref_state;
    delete ref_control;
    delete controller;
  }

  void run() {

    switch(state){
      case RunState::Idle:

      break;
      case RunState::Init:
        if(turnToAngle(ref_state->get(0))){
          clock_gettime(CLOCK_REALTIME, &time_since_last);
          run_counter = 0;
          state = RunState::Running;
        }
      break;
      case RunState::Running:
        clock_gettime(CLOCK_REALTIME, &time_now);
        unsigned long elapsed_t = time_now.tv_nsec - time_since_last.tv_nsec;

        if(elapsed_t > nano_dt){
          if(!computeControl()){
            cmd_u << 0.0, 0.0;
            state = RunState::Idle;
          }
          run_counter++;
          clock_gettime(CLOCK_REALTIME, &time_since_last);
        }
      break;
    }

    geometry_msgs::Twist msg;
    msg.linear.x = cmd_u[0];
    msg.angular.z = cmd_u[1];
    movement_pub.publish(msg);
  }

  bool turnToAngle(Eigen::Vector3d vec){
    double angle_diff = mpc::Tools::subtractAngles(vec[2], pose.getTheta());
    ROS_INFO("Turning to first point, angle difference = %.2f", angle_diff);
    if(angle_diff < -0.1){
      cmd_u << 0.0, -1.5;
      return false;
    }
    else if(angle_diff > 0.1){
      cmd_u << 0.0, 1.5;
      return false;
    }
    return true;
  }

  bool computeControl(){
    mpc::StateVectors* horizon_x;
    mpc::ControlVectors* horizon_u;
    bool failed = false;

    int horizon_length = 15;
    if(this->trajectory.size() <= run_counter + horizon_length)
      horizon_length = this->trajectory.size() - run_counter - 1;

    if(horizon_length < 2)
      return false;

    try{
      horizon_x = ref_state->getHorizon(run_counter, horizon_length);
    } catch(std::exception& e){
      std::cout << "State horizon exception." << std::endl;
      return false;
    }

    try{
      horizon_u = ref_control->getHorizon(run_counter, horizon_length);
    } catch(std::exception& e){
      std::cout << "Control horizon exception." << std::endl;
      return false;
    }

    ROS_INFO("Horizon x size = %d, y size = %d", (int)horizon_x->size(), (int)horizon_u->size());
    try{
      cmd_u = controller->computeControl(horizon_x, horizon_u, pose.getPoseVector(), dt);
    } catch(std::exception& e){
      std::cout << "Compute Control exception." << std::endl;
      return false;
    }
    
    std::cout << "Command u = " << cmd_u << std::endl;
    delete horizon_x;
    delete horizon_u;
    return true;
  }  

  bool ok() { return nh.ok(); }
};

int main(int argc, char *argv[])  {

  // Init the connection with the ROS system.
  ros::init(argc, argv, "control_controller_node");

  // Create a class instance.
  Robot robot(1.65, 2.5, 1.5);
  
  // Start the ROS main loop and check for incoming sensor messages.
  ros::Rate loop_rate(30);
  while (robot.ok()) {

    // Call our own procedures during the idle time of the ROS loop.
    robot.run();

    // Sleep for the remaining idle time.
    ros::spinOnce();
    loop_rate.sleep();
  }  
  return 0;
}