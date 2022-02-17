#include "Headers/Common.hpp"
class PointCloudProcessor {
    
    // Given a PointCloud, process it to push to the buffer

    public:
        Points points;
        PointCloudProcessor(const PointCloud_msg&);

    private:
        Points msg2points(const PointCloud_msg&);

        // Velodyne specific
            Points velodynemsg2points(const PointCloud_msg&);
            double get_begin_time(const pcl::PointCloud<velodyne_ros::Point>&);
        
        // HESAI specific
            Points hesaimsg2points(const PointCloud_msg&);
            double get_begin_time(const pcl::PointCloud<hesai_ros::Point>&);
        
        // Ouster specific
            Points oustermsg2points(const PointCloud_msg&);
            double get_begin_time(const pcl::PointCloud<ouster_ros::Point>&);

        // Custom specific
            Points custommsg2points(const PointCloud_msg&);
            double get_begin_time(const pcl::PointCloud<custom::Point>&);
        
        template <typename PointType> Points to_points(const typename pcl::PointCloud<PointType>&);

        Points downsample(const Points&);
        static bool time_sort(const Point&, const Point&);
        void sort_points(Points&);
};
