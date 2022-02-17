#ifndef __OBJECTS_H__
#define __OBJECTS_H__
#include "Headers/Common.hpp"
#include "Headers/Utils.hpp"
#include "Headers/Objects.hpp"
#include "Headers/Publishers.hpp"
#include "Headers/PointClouds.hpp"
#include "Headers/Accumulator.hpp"
#include "Headers/Compensator.hpp"
#include "Headers/Localizator.hpp"
#include "Headers/Mapper.hpp"
#endif

extern struct Params Config;

// class Accumulator
    // public:
        // Add content to buffer
        void Accumulator::add(State cnt, double time) {
            if (time > 0) cnt.time = time;
            this->push(cnt);
        }

        void Accumulator::add(IMU cnt, double time) {
            if (time > 0) cnt.time = time;
            this->push(cnt);
        }

        void Accumulator::add(Point cnt, double time) {
            if (time > 0) cnt.time = time;
            this->push(cnt);
        }

        // 拷贝点云
        void Accumulator::add(Points points) {
            for (Point p : points) this->push(p);
        }

        // Receive from topics
        // 激光雷达回调
        void Accumulator::receive_lidar(const PointCloud_msg& msg) {
            // 创建临时对象PointCloudProcessor，将点云转换为自定义类型Point，并且完成降采样
            PointCloudProcessor processed(msg);
            // 拷贝点云到Accumulator::BUFFER_L
            this->add(processed.points);
        }

        // msg类型转换，拷贝imu数据到Accumulator::BUFFER_I
        void Accumulator::receive_imu(const IMU_msg& msg) {
            IMU imu(msg);
            this->add(imu);
        }

        // Empty buffers
        void Accumulator::clear_buffers() {
            this->BUFFER_L.clear();
            this->BUFFER_I.clear();
        }

        void Accumulator::clear_buffers(TimeType t) {
            this->BUFFER_L.clear(t);
            this->BUFFER_I.clear(t);
        }

        void Accumulator::clear_lidar(TimeType t) {
            this->BUFFER_L.clear(t);
        }

        /////////////////////////////////

        State Accumulator::get_prev_state(double t) {
            return this->get_prev(this->BUFFER_X, t);
        }

        IMU Accumulator::get_next_imu(double t) {
            return this->get_next(this->BUFFER_I, t);
        }

        States Accumulator::get_states(double t1, double t2) {
            return this->get(this->BUFFER_X, t1, t2);
        }

        Points Accumulator::get_points(double t1, double t2) {
            return this->get(this->BUFFER_L, t1, t2);
        }

        IMUs Accumulator::get_imus(double t1, double t2) {
            return this->get(this->BUFFER_I, t1, t2);
        }

        //////////////////////////

        bool Accumulator::ready() {
            // Only check it once
            if (this->is_ready) return true;
            
            // Ready if there's enough IMUs to fit the delay
            // 如果有足够的imu数据，is_ready设置为true，并且设置initial_time
            if (this->enough_imus()) {
                this->set_initial_time();
                return this->is_ready = true;
            }

            return this->is_ready = false;
        }

        bool Accumulator::ended(double t) {
            if (not this->ready()) return false;
            if (t - initial_time < 3) return false;
            return this->get_imus(t - 3., t).size() < 2;
        }

        ros::Rate Accumulator::refine_delta(const HeuristicParams& heuristic, double t) {
            assert(("There has to be exactly one more delta value than time delimiters", heuristic.times.size() + 1 == heuristic.deltas.size()));
            this->delta = this->interpret_heuristic(heuristic, t);
            return ros::Rate((int) 1./this->delta);
        }

    // private:

        void Accumulator::push(const State& state) { this->BUFFER_X.push(state); }
        void Accumulator::push(const IMU& imu) { this->BUFFER_I.push(imu); }
        void Accumulator::push(const Point& point) { this->BUFFER_L.push(point); }

        bool Accumulator::enough_imus() {
            return this->BUFFER_I.size() > Config.real_time_delay*Config.imu_rate;
        }

        void Accumulator::set_initial_time() {
            if (this->BUFFER_I.size() < 1) return;
            double latest_imu_time = this->BUFFER_I.front().time;

            this->initial_time = latest_imu_time - Config.real_time_delay;
        }

        double Accumulator::interpret_heuristic(const HeuristicParams& heuristic, double t) {
            // If is after last time
            if (heuristic.times.empty()) return heuristic.deltas.back();
            if (t - this->initial_time >= heuristic.times.back()) return heuristic.deltas.back();
            
            // If we have to find it in the list
            for (int k = 0; k < heuristic.times.size(); ++k)
                if (t - this->initial_time < heuristic.times[k])
                    return heuristic.deltas[k];

            return heuristic.deltas.back();
        }
