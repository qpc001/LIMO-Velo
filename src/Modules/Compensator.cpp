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

// class Compensator
    // public:
        Points Compensator::compensate(double t1, double t2) {
            // Call Accumulator
            Accumulator& accum = Accumulator::getInstance();

            // Points from t1 to t2
            Points points = accum.get_points(t1, t2);   // 获取[t1,t2]时间段内的点
            if (points.empty()) return Points();

            // (Integrated) States surrounding t1 and t2
            // 得到了before t1 and after t2的状态，时间区间为[ < t1, > t2]
            States path_taken = this->path(t1, t2);
            assert (path_taken.size() >= 2);
            ROS_WARN("Compensator::compensate path States size[%zu]", path_taken.size());

            // Compensated points given a path
            // 得到t2时刻的状态
            State Xt2 = this->get_t2(path_taken, t2);

            //
            return this->compensate(path_taken, Xt2, points);
        }

        /**
         * @brief Compensator::path
         * @param t1
         * @param t2
         * @return 得到了before t1 and after t2的状态，时间区间为[ < t1, > t2]
         */
        States Compensator::path(double t1, double t2) {
            // Call Accumulator
            Accumulator& accum = Accumulator::getInstance();

            // Get states just before t1 to t2
            // 获取[t1,t2]时间段内的所有状态
            States states = accum.get_states(t1, t2);
            ROS_WARN("Compensator::path States size[%zu]", states.size());
            // 然后在最前面加上t1时刻的前一个状态
            states.push_front(accum.get_prev_state(t1));

            // Get imus from first state to just after t2
            // 获取所有states[]内时间段的imu数据
            IMUs imus = accum.get_imus(states.front().time, t2);
            imus.push_back(accum.get_next_imu(t2));

            // 得到了before t1 and after t2的状态，时间区间为[ < t1, > t2]
            return this->upsample(states, imus);
        }

    // private:
        /**
         * @brief Compensator::get_t2
         * 得到t2时刻的状态(使用1帧imu推算得到)
         * @param states
         * @param t2
         * @return
         */
        State Compensator::get_t2(const States& states, double t2) {
            // 首先指向最后一个状态
            int s = states.size() - 1;
            assert (states.front().time <= t2);
            // 如果当前指向的状态时间戳 > t2, 就指向前一个状态
            while (t2 < states[s].time) --s;
            
            // 此时 Xt2 时间戳 <= t2
            State Xt2 = states[s];
            // 构造了一个IMU对象
            // 进行状态推算，推算到t2时刻
            Xt2 += IMU (Xt2.a, Xt2.w, t2);
            return Xt2;
        }

        /*
            @Input:
                states: before t1 and to t2
                imus: before t1 and after t2
            
            @Output:
                upsampled_states (size := imus.size): before t1 and after t2
        */
        States Compensator::upsample(const States& states, const IMUs& imus) {
            // 确保第一个imu时间< 第一个状态时间  ，并且 ，  最后一个状态时间 < 最后一个 imu时间
            assert (imus.front().time <= states.front().time and states.back().time <= imus.back().time);

            int s, u;
            s = u = 0;

            States upsampled_states;
            State int_state = states[s];    //

            // IMUs between two states
            while (s < states.size() - 1) { // 遍历： 如果不是最后一个状态
                // 先保存初始状态到upsampled_states[]
                upsampled_states.push_back(states[s]);

                // 遍历imu， 取下一个状态之前的imu，进行推算
                while (u < imus.size() and imus[u].time < states[s+1].time) {
                    // 这里重载了 + 运算符，实际操作是使用 imu数据来对 int_state 状态进行推算
                    int_state += imus[u++];
                    // 每推算一次imu，就将状态保存到upsampled_states
                    upsampled_states.push_back(int_state);
                }

                // 推算完之后，int_state 直接取下一个状态
                int_state = states[s++];
            }

            if (u >= imus.size()) u = imus.size() - 1;  // 避免u越界
            // 将最后一个状态保存到upsampled_states[]
            upsampled_states.push_back(states.back());
            int_state = states.back();  // int_state取最后一个状态

            // IMUs after last state
            // 如果最后一个imu的时间戳 > 最后一个状态
            while (int_state.time < imus.back().time) {
                // 为啥是循环， 这里有可能有越界的风险
                int_state += imus[u++]; // 状态推算
                upsampled_states.push_back(int_state);  // 保存状态
            }

            // 到这里，得到了before t1 and after t2的状态，时间区间为[ < t1, > t2]
            return upsampled_states;
        }

        Points Compensator::downsample(const Points& points) {
            return this->voxelgrid_downsample(points);
            // return this->onion_downsample(points);
        }

        /*
            @Input:
                states: path the car has taken (pre and post included)
                points: stamped points during path
            @Output:
                compensated_points: compensated points ready to be transported by Xt2
            
            @Pseudocode:
                for each state:
                    for each point between state and next_state:
                        compensate point matching its time via integrating state's last IMU
        */

        /**
         * @brief Compensator::compensate
         * 将每一个点投影到 t2时刻下的lidar坐标系
         * @param states
         * @param Xt2
         * @param points
         * @return
         */
        Points Compensator::compensate(const States& states, const State& Xt2, const Points& points) {
            // States have to surround points
            // 确保 第一个状态时间戳 <= 第一个点时间戳， 最后一个状态时间戳 >= 最后一个点时间戳
            assert (not states.empty() and states.front().time <= points.front().time and  points.back().time <= states.back().time);

            Points t2_inv_ps;
            int p = 0;

            // 遍历状态
            for (int s = 0; s < states.size() - 1; ++s) {
                // 遍历点集： 如果当前状态时间 < 当前点时间 ， 并且 当前点时间 < 下一个状态时间，即可进行畸变矫正
                while (p < points.size() and states[s].time <= points[p].time and points[p].time <= states[s+1].time) {                    
                    // Integrate to point time
                    // 取当前状态
                    State Xtp = states[s];
                    // 推算到该points点的时间戳
                    Xtp += IMU (states[s].a, states[s].w, points[p].time);

                    // Transport to X_t2^-1 frame
                    // Xtp: 该points点的时间戳对应的pose（imu）
                    // Xtp.I_Rt_L: 激光雷达到imu的变换
                    // global_p: 该points点的时间戳在世界坐标系的坐标 (推算得到)
                    Point global_p = Xtp * Xtp.I_Rt_L() * points[p];
                    // Xt2: t2时刻下的pose(imu)
                    // Xt2.inv() * global_p: 表示将世界坐标系的点投影到 t2时刻下的imu坐标系
                    // t2_inv_p: 投影到 t2时刻下的lidar坐标系的点
                    // 完成畸变矫正
                    Point t2_inv_p = Xt2.I_Rt_L().inv() * Xt2.inv() * global_p;
                    t2_inv_ps.push_back(t2_inv_p);

                    ++p;
                }
            }

            return t2_inv_ps;
        }

        Points Compensator::voxelgrid_downsample(const Points& points) {
            // Create a PointCloud pointer
            pcl::PointCloud<full_info::Point>::Ptr pcl_ptr(new pcl::PointCloud<full_info::Point>());
            Processor::fill(*pcl_ptr, points);

            // Downsample using a VoxelGrid
            pcl::PointCloud<full_info::Point> ds_pcl;
            pcl::VoxelGrid<full_info::Point> filter;
            filter.setInputCloud(pcl_ptr);
            filter.setLeafSize(0.5, 0.5, 0.5);
            filter.filter(ds_pcl);
            
            Points ds_points;
            for (auto p : ds_pcl.points) ds_points.push_back(Point (p));
            return ds_points;
        }

        Points Compensator::onion_downsample(const Points& points) {
            Points ds_points;

            for (int i = 0; i < points.size(); ++i) {
                const Point& p = points[i];
                if (0 < p.range and p.range < 4 and (256/Config.ds_rate <= 1 or i%(256/Config.ds_rate) == 0)) ds_points.push_back(p);
                else if (4 < p.range and p.range < 6 and (64/Config.ds_rate <= 1 or i%(64/Config.ds_rate) == 0)) ds_points.push_back(p);
                else if (6 < p.range and p.range < 9 and (32/Config.ds_rate <= 1 or i%(32/Config.ds_rate) == 0)) ds_points.push_back(p);
                else if (9 < p.range and p.range < 12 and (16/Config.ds_rate <= 1 or i%(16/Config.ds_rate) == 0)) ds_points.push_back(p);
                else if (12 < p.range and p.range < 22 and (8/Config.ds_rate <= 1 or i%(8/Config.ds_rate) == 0)) ds_points.push_back(p);
                else if (22 < p.range and p.range < 30 and (4/Config.ds_rate <= 1 or i%(4/Config.ds_rate) == 0)) ds_points.push_back(p);
                else if (30 < p.range and p.range < 50 and (2/Config.ds_rate <= 1 or i%(2/Config.ds_rate) == 0)) ds_points.push_back(p);
                else if (p.range > 50) ds_points.push_back(p);
            }

            return ds_points;
        }
