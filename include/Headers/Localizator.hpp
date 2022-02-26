#ifndef SLAM_LIBRARIES_H
#define SLAM_LIBRARIES_H
// SLAM libraries
#include "use-ikfom.hpp"
#include "ikd_Tree.h"
#endif

class Localizator {
    
    // Variables

    public:
        Points points2match;
        double last_time_integrated = -1;
        double last_time_updated = -1;

    private:
        esekfom::esekf<state_ikfom, 12, input_ikfom> IKFoM_KF;
        State initial_state_;
    // Methods

    public:
        Localizator();
        void update(const Points&, double time);
        void calculate_H(const state_ikfom&, const Matches&, Eigen::MatrixXd& H, Eigen::VectorXd& h);
        
        void propagate_to(double t1, double t);
        State initial_state();
        State latest_state();

    private:
        void init_IKFoM();
        void init_IKFoM_state();
        void init_IKFoM_orientation(const IMU& imu);
        void IKFoM_update(const Points&);
        
        void propagate(const IMU& imu);
        const state_ikfom& get_x() const;

    // Singleton pattern

    public:
        static Localizator& getInstance();
    
    private:
        // Delete copy/move so extra instances can't be created/moved.
        Localizator(const Localizator&) = delete;
        Localizator& operator=(const Localizator&) = delete;
        Localizator(Localizator&&) = delete;
        Localizator& operator=(Localizator&&) = delete;

};