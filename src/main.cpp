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

Params Config;

int main(int argc, char** argv) {
    ros::init(argc, argv, "limovelo");
    ros::NodeHandle nh;

    // Get YAML parameters
    nh.param<bool>("mapping_online", Config.mapping_online, true);
    nh.param<bool>("real_time", Config.real_time, true);
    nh.param<bool>("estimate_extrinsics", Config.estimate_extrinsics, false);
    nh.param<bool>("print_extrinsics", Config.print_extrinsics, false);
    nh.param<int>("ds_rate", Config.ds_rate, 4);
    nh.param<int>("MAX_NUM_ITERS", Config.MAX_NUM_ITERS, 3);
    nh.param<std::vector<double>>("LIMITS", Config.LIMITS, std::vector<double> (23, 0.001));
    nh.param<int>("NUM_MATCH_POINTS", Config.NUM_MATCH_POINTS, 5);
    nh.param<int>("MAX_POINTS2MATCH", Config.MAX_POINTS2MATCH, 10);
    nh.param<double>("MAX_DIST_PLANE", Config.MAX_DIST_PLANE, 2.0);
    nh.param<float>("PLANES_THRESHOLD", Config.PLANES_THRESHOLD, 0.1f);
    nh.param<float>("PLANES_CHOOSE_CONSTANT", Config.PLANES_CHOOSE_CONSTANT, 9.0f);
    nh.param<std::string>("LiDAR_type", Config.LiDAR_type, "unknown");
    nh.param<double>("LiDAR_noise", Config.LiDAR_noise, 0.001);
    nh.param<double>("min_dist", Config.min_dist, 3.);
    nh.param<double>("imu_rate", Config.imu_rate, 400);
    nh.param<double>("degeneracy_threshold", Config.degeneracy_threshold, 5.);
    nh.param<bool>("print_degeneracy_values", Config.print_degeneracy_values, false);
    nh.param<double>("full_rotation_time", Config.full_rotation_time, 0.1);
    nh.param<double>("empty_lidar_time", Config.empty_lidar_time, 20.);
    nh.param<double>("real_time_delay", Config.real_time_delay, 1.);
    nh.param<double>("covariance_gyroscope", Config.cov_gyro, 1e-4);
    nh.param<double>("covariance_acceleration", Config.cov_acc, 1e-2);
    nh.param<double>("covariance_bias_gyroscope", Config.cov_bias_gyro, 1e-5);
    nh.param<double>("covariance_bias_acceleration", Config.cov_bias_acc, 1e-4);
    nh.param<double>("wx_MULTIPLIER", Config.wx_MULTIPLIER, 1);
    nh.param<double>("wy_MULTIPLIER", Config.wy_MULTIPLIER, 1);
    nh.param<double>("wz_MULTIPLIER", Config.wz_MULTIPLIER, 1);
    nh.param<std::string>("points_topic", Config.points_topic, "/velodyne_points");
    nh.param<std::string>("imus_topic", Config.imus_topic, "/vectornav/IMU");
    nh.param<std::vector<double>>("/Heuristic/times", Config.Heuristic.times, {1.});
    nh.param<std::vector<double>>("/Heuristic/deltas", Config.Heuristic.deltas, {0.1, 0.01});
    nh.param<std::vector<float>>("initial_gravity", Config.initial_gravity, {0.0, 0.0, -9.807});
    nh.param<std::vector<float>>("I_Translation_L", Config.I_Translation_L, std::vector<float> (3, 0.));
    nh.param<std::vector<float>>("I_Rotation_L", Config.I_Rotation_L, std::vector<float> (9, 0.));

    // Objects
    Publishers publish(nh); // 创建一系列publisher
    // 创建几个模块，单例模式
    Accumulator& accum = Accumulator::getInstance();
    Compensator comp = Compensator ();
    Mapper& map = Mapper::getInstance();
    Localizator& KF = Localizator::getInstance();

    // Subscribers
    // 订阅激光、IMU
    ros::Subscriber lidar_sub = nh.subscribe(
        Config.points_topic, 1000,
        &Accumulator::receive_lidar, &accum
    );

    ros::Subscriber imu_sub = nh.subscribe(
        Config.imus_topic, 1000,
        &Accumulator::receive_imu, &accum
    );

    ros::Rate rate(10);

    // If not real time, define an artificial clock
    double clock = -1;
    double t2 = -1;

    while (ros::ok()) {
        // 如果有足够的imu数据，is_ready设置为true，会进入循环
        while (accum.ready()) {
            // 是否需要实时？
            if (Config.real_time) {
                // Should be t2 = ros::Time::now() - delay
                double latest_imu_time = accum.BUFFER_I.front().time;
                t2 = latest_imu_time - Config.real_time_delay; 
            } else {
                // 刚开始时，clock = -1
                if (clock < 0) clock = accum.initial_time;  // 给clock赋值
                t2 = clock - Config.real_time_delay; 
                clock += accum.delta;
            }

            // Refine delta if need to be
            rate = accum.refine_delta(Config.Heuristic, t2);
            
            // Define t1
            double t1 = std::max(KF.last_time_updated, t2 - accum.delta);   // t1: 上一次更新的时间?
            if (KF.last_time_updated < 0) t1 = t2 - accum.delta;    // 如果还没有优化过，即第一次，直接赋值t1= t2 - accum.delta

            // Integrate from t1 to t2
            KF.propagate_to(t2);    // imu传播

            // Field of view too small (relevant for real-time)
            if (t2 - t1 < accum.delta - 1e-6) break;    // 如果t2 - t1太小了，时间太短，激光数据不够

            // 建图 或者 定位
            if (Config.mapping_online or (not Config.mapping_online and map.exists())) {
                // Compensated pointcloud given a path
                // 去畸变， 将每一个点投影到 t2时刻下的lidar坐标系
                Points compensated = comp.compensate(t1, t2);
                // 0.5, 0.5, 0.5的voxel降采样
                Points ds_compensated = comp.downsample(compensated);
                if (ds_compensated.size() < Config.MAX_POINTS2MATCH) break; 

                // Localize points in map
                // 使用去畸变后的点云进行状态更新
                KF.update(ds_compensated, t2);
                State Xt2 = KF.latest_state();
                accum.add(Xt2, t2);
                publish.state(Xt2, false);
                publish.tf(Xt2);

                // Publish compensated
                Points global_compensated = Xt2 * Xt2.I_Rt_L() * ds_compensated;
                publish.pointcloud(global_compensated, true);

                // Map at the same time (online)
                if (Config.mapping_online) {
                    map.add(global_compensated, t2, true);
                    publish.pointcloud(global_compensated, false);
                    
                    if (Config.print_extrinsics) publish.extrinsics(Xt2);
                }
            }
            
            // Add updated points to map (offline)
            if (not Config.mapping_online and map.hasToMap(t2)) {
                State Xt2 = KF.latest_state();
                Points full_compensated = comp.compensate(t2 - Config.full_rotation_time, t2);
                Points global_full_compensated = Xt2 * Xt2.I_Rt_L() * full_compensated;
                Points global_full_ds_compensated = comp.downsample(global_full_compensated);

                if (global_full_ds_compensated.size() < Config.MAX_POINTS2MATCH) break; 
                if (Config.print_extrinsics) publish.extrinsics(Xt2);

                map.add(global_full_ds_compensated, t2, true);
                publish.pointcloud(global_full_ds_compensated, false);
            }

            // Empty too old LiDAR points
            accum.clear_lidar(t2 - Config.empty_lidar_time);

            // Trick to call break in the middle of the program
            break;
        }

        if (accum.ended(t2) and not Config.real_time) break;

        ros::spinOnce();
        rate.sleep();
    }

    ROS_ERROR("Rosbag ending detected.");

    return 0;
}
